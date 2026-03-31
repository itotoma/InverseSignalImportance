# ISI.py
"""
Inverse Signal Importance (ISI) — standalone modeling module.

This module provides tools for analyzing and predicting time-series data
using a hybrid EM algorithm that combines Ridge regression with Kalman smoothing.

Key features:
- Continuous-time Kalman filter and RTS smoother
- Ridge regression for response function estimation
- EM algorithm for iterative parameter optimization
- Multi-group (cross-validation-style) training support

Example
-------
    import ISI
    em = ISI.EMAlgorithm(window_size=3, process_noise_var=10.0)
    x_smooth, P_smooth, coef, rss, rss_history = em.fit(
        X_list=[X1, X2], Y_list=[Y1, Y2], times=time_array, n_signals=4
    )
    test_rss, predictions = em.test(
        X_list=[X3], Y_list=[Y3], x_smooth=x_smooth, coef=coef
    )
"""


import numpy as np
import pandas as pd
from typing import Sequence, Tuple, Optional
from scipy.integrate import solve_ivp
from scipy import interpolate
from sklearn.linear_model import Ridge


class EMAlgorithm:

    """
    Hybrid EM algorithm combining Ridge regression and Kalman smoothing
    for nonlinear time-varying signal importance estimation.

    At each EM iteration:
    - step 1: estimate latent signal importance states via Kalman smoothing
    - step 2: update response function coefficients via Ridge regression

    Parameters
    ----------
    window_size : int
        Number of consecutive days in each sliding input window.
        Controls the length of the response function per signal.
    process_noise_var : float
        Variance of the process noise (Q = I * process_noise_var).
        Larger values allow signal importance to vary more rapidly over time.
    observe_noise_var : float
        Diagonal variance of the observation noise matrix R.
        Should be set based on the variance of dGSI across data groups.
    observe_noise_cov : float
        Off-diagonal covariance of R, shared across data groups.
        Should be set based on the covariance of dGSI between groups.
    alpha : float
        Ridge regularization strength for response function estimation.
    """

    def __init__(
        self,
        window_size: int,
        process_noise_var: float,
        observe_noise_var: float = 1.5,
        observe_noise_cov: float = 0.35,
        alpha: float = 1.0
    ):
        self.window_size = window_size
        self.alpha = alpha
        self.process_noise_var = process_noise_var
        self.observe_noise_var = observe_noise_var
        self.observe_noise_cov = observe_noise_cov

        # Number of Kalman smoothing passes per EM iteration.
        # Multiple passes stabilize state estimates before the Ridge update.
        self.kalman_repeat = 3
        
        self.F: Optional[np.ndarray] = None  # State transition matrix
        self.W: Optional[np.ndarray] = None  # Process noise covariance
        self.R: Optional[np.ndarray] = None  # Observation noise covariance

    def fit(
        self,
        X_list: pd.DataFrame, # Cross validation list, given as stack form
        Y_list: pd.DataFrame, # Cross validation list, given as stack form
        times: np.ndarray,
        n_signals: int,
        max_iter: int = 1000,
        tol: float = 1e-7
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray]:
        
        """
        Run the EM algorithm until convergence.

        Parameters
        ----------
        X_list : list of np.ndarray, each shape (T, window_size * n_signals)
            Sliding-window input matrices, one per data group.
        Y_list : list of np.ndarray, each shape (T,)
            Observed output (dGSI) time series, one per data group.
        times : np.ndarray, shape (T,)
            Observation time points (in raw day units).
        n_signals : int
            Number of environmental signals (equals the latent state dimension).
        max_iter : int
            Maximum number of EM iterations.
        tol : float
            Convergence threshold on RSS improvement.

        Returns
        -------
        x_smooth : np.ndarray, shape (T, n_signals)
            Smoothed signal importance estimates over time.
        P_smooth : np.ndarray, shape (T, n_signals, n_signals)
            Smoothed state covariance matrices.
        coef : np.ndarray, shape (n_signals * window_size,)
            Estimated response function coefficients.
        rss : float
            Final training RSS.
        rss_history : np.ndarray
            RSS at each EM iteration (truncated at convergence).
        """

        # validate input data structure 
        state_dim = n_signals
        sample_num = len(Y_list) 
        if len(Y_list) != len(X_list):
            raise ValueError(f"X and Y must have same number of samples, got {X_stack.shape} and {Y_stack.shape}")
        if len(times) < 2:
            raise ValueError("times must have at least 2 elements")
        if state_dim <= 0:
            raise ValueError("n_signals must be positive")
        if max_iter <= 0:
            raise ValueError("max_iter must be positive")
        if tol <= 0:
            raise ValueError("tol must be positive")

    
        # Stack groups for Ridge regression
        X_stack = np.vstack(X_list)
        Y_stack = np.array(Y_list).reshape(-1,1)
        Y = np.array(Y_list).T
        dt = np.diff(times)/46
        T = len(times)
        
 
        # Initialize
        self.F = np.eye(state_dim)
        self.W = np.eye(state_dim) * self.process_noise_var
        # R: observation noise covariance across groups
        # Diagonal = observe_noise_var, off-diagonal = observe_noise_cov
        self.R = np.ones((sample_num)) * self.observe_noise_cov
        self.R = self.R - np.eye(sample_num)*self.observe_noise_cov + np.eye(sample_num)*self.observe_noise_var
        x0 = np.ones(state_dim)
        P0 = np.eye(state_dim) * 100

        # EM loop
        X_update = X_stack.copy()
        prev_rss = np.inf
        rss_history = np.zeros(max_iter)
        
        x_filt = np.zeros((T, state_dim))
        P_filt = np.zeros((T, state_dim, state_dim))
        x_smooth = np.zeros((T, state_dim))
        P_smooth = np.zeros((T, state_dim, state_dim))
        
        for k in range(max_iter):
            
            # ---M-step: update response function coefficients ---
            if self.window_size == 1:
                H_list = [X for X in X_list]
                coef = np.ones(state_dim)
            else:
                coef = self._fit_ridge(X_update, Y_stack, self.alpha)
                transA = np.zeros((state_dim * self.window_size, state_dim))
                for i in range(1, state_dim + 1):
                    start = self.window_size * (i - 1)
                    end = self.window_size * i
                    transA[start:end, i-1] = coef[start:end]
                H_list = [X @ transA for X in X_list]
            H = np.stack([np.vstack([Hl[h] for Hl in H_list]) for h in range(T)])
            
            # --- E-step: Kalman smoothing (multiple passes for stability) ---
            for _ in range(self.kalman_repeat):
                x_filt, P_filt, x_smooth, P_smooth = self._kalman_smooth(
                    H, Y, x0, P0, dt
                )
                x0 = x_smooth[-1]
                P0 = P_filt[-1]
            
            # --- Compute training RSS ---
            train_pred_list = np.array([[Hl[t] @ x_smooth[t] for t in range(T)] for Hl in H_list])
            rss = np.mean(np.sum((Y_list - train_pred_list)**2, axis=1))
            rss_history[k] = rss
            
            # --- Check convergence ---
            diff = prev_rss - rss
            prev_rss = rss
            if diff < tol:
                print(f"rss: {rss}")
                rss_history = rss_history[:k+1]
                break

            # --- Update X for next Ridge step ---
            ExtendedZ = np.repeat(x_smooth.T, self.window_size, axis=0)
            X_update = np.vstack([X * ExtendedZ.T for X in X_list])
    
        return x_smooth, P_smooth, coef, rss, rss_history
    
    
    def _fit_ridge(self, X: np.ndarray, Y: np.ndarray, alpha: float) -> np.ndarray:
        """Fit Ridge regression and return coefficient matrix (supports multi-output)."""
        if X.shape[0] != Y.shape[0]:
            raise ValueError(f"X and Y must have same number of samples, got {X.shape} and {Y.shape}")
        if alpha < 0:
            raise ValueError("alpha must be non-negative")
        if np.any(np.isnan(X)) or np.any(np.isnan(Y)):
            raise ValueError("Input data contains NaN values")
        if np.any(np.isinf(X)) or np.any(np.isinf(Y)):
            raise ValueError("Input data contains infinite values")
        model = Ridge(alpha=alpha)
        model.fit(X, Y)
        return model.coef_[0]
    
    
    def _kalman_filter(
        self,
        x: np.ndarray,
        P: np.ndarray,
        y: np.ndarray,
        F: np.ndarray,
        W: np.ndarray,
        H: np.ndarray,
        R: np.ndarray,
        delta: float
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Continuous-time Kalman filter: predict and update."""
        # Prediction
        sol = solve_ivp(
            lambda t, P_flat: (self.F @ P_flat.reshape(self.F.shape) @ self.F.T + self.W).flatten(),
            (0, delta), P.flatten(), t_eval=[delta]
        )
        x_pred = self.F @ x
        P_pred = sol.y[:, -1].reshape(P.shape)
        
        # Update
        S = H @ P_pred @ H.T + R
        y_pred = H @ x_pred
        
        K = (np.linalg.solve(S.T, H @ P_pred.T)).T 
        x_upd = x_pred + K @ (y - y_pred)
        P_upd = P_pred - K @ H @ P_pred
        
        return x_upd, P_upd, x_pred, P_pred, y_pred, S
    
    
    def _rts_smoother(
        self,
        x_filt: np.ndarray,
        P_filt: np.ndarray,
        F: np.ndarray,
        W: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Rauch–Tung–Striebel smoother."""
        T, n = x_filt.shape
        x_smooth = x_filt.copy()
        P_smooth = P_filt.copy()
        
        # 逆行列の計算を最適化
        for t in range(T - 2, -1, -1):
            P_pred = F @ P_smooth[t] @ F.T + W
            
            K = np.linalg.solve(P_pred.T, F @ P_smooth[t].T).T
            x_smooth[t] += K @ (x_smooth[t + 1] - F @ x_smooth[t])
            P_smooth[t] += K @ (P_smooth[t + 1] - P_pred) @ K.T
            
        return x_smooth, P_smooth
    

    def _kalman_smooth(
        self,
        H: np.ndarray,
        Y: np.ndarray,
        x0: np.ndarray,
        P0: np.ndarray,
        dt: np.ndarray
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        T = Y.shape[0]
        if T != len(dt) + 1:
            raise ValueError(f"Length of Y ({T}) must be equal to length of dt + 1 ({len(dt) + 1})")
        if x0.shape[0] != self.F.shape[0]:
            raise ValueError(f"Dimension of x0 ({x0.shape[0]}) must match F ({self.F.shape[0]})")
        if P0.shape != (self.F.shape[0], self.F.shape[0]):
            raise ValueError(f"Shape of P0 {P0.shape} must match F {self.F.shape}")
        x_filt = np.zeros((T, self.F.shape[0]))
        P_filt = np.zeros((T, self.F.shape[0], self.F.shape[1]))
        for t in range(T):
            if t == 0:
                x_filt[t], P_filt[t],_,_,_,_  = self._kalman_filter(
                    x0, P0, Y[t,], self.F, self.W, H[t], self.R, 1
                )
            else:
                x_filt[t], P_filt[t],_,_,_,_  = self._kalman_filter(
                    x_filt[t-1], P_filt[t-1], Y[t,], self.F, self.W, H[t], self.R, dt[t-1]
                )
      
        x_smooth, P_smooth = self._rts_smoother(x_filt, P_filt, self.F, self.W)
        return x_filt, P_filt, x_smooth, P_smooth

    def test(
        self,
        X_list: pd.DataFrame,
        Y_list: pd.DataFrame,
        x_smooth: np.ndarray,
        coef: np.ndarray
    ) -> [np.ndarray, np.ndarray]:
        """
        Evaluate the model on held-out data using learned states and coefficients.

        Parameters
        ----------
        X_list : list of np.ndarray, each shape (T, window_size * n_signals)
            Input matrices for the test group(s).
        Y_list : list of np.ndarray, each shape (T,)
            Observed output for the test group(s).
        x_smooth : np.ndarray, shape (T, n_signals)
            Smoothed states from training (shared across groups).
        coef : np.ndarray, shape (n_signals * window_size,)
            Response function coefficients from training.

        Returns
        -------
        rss_test : float
            Mean RSS across test groups.
        test_pred : np.ndarray, shape (n_test_groups, T)
            Predicted outputs for each test group.
        """
        if coef.shape[0] != x_smooth.shape[1] * self.window_size:
            raise ValueError(f"coef shape mismatch: expected {x_smooth.shape[1] * self.window_size}, got {coef.shape[0]}")

        
        transA = np.zeros((x_smooth.shape[1] * self.window_size, x_smooth.shape[1]))
        for i in range(1, x_smooth.shape[1] + 1):
            start = self.window_size * (i - 1)
            end = self.window_size * i
            transA[start:end, i-1] = coef[start:end]
        
        # 予測値の計算
        H_list = [X @ transA for X in X_list]
        test_pred_list = np.array([[Hl[t] @ x_smooth[t] for t in range(len(x_smooth))] for Hl in H_list])
        rss_test = np.mean(np.sum((Y_list - test_pred_list)**2, axis=1))
        
        return rss_test, test_pred_list


