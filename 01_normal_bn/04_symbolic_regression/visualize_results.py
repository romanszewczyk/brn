"""
Visualization script for PySR symbolic regression results.
"""

import numpy as np
import scipy.io
import matplotlib.pyplot as plt

# Load data
mat_data = scipy.io.loadmat('optimal_bn_500times.mat')
n = mat_data['n'].flatten()
bn_all = mat_data['bn_all']
bn_avg = np.mean(bn_all, axis=0)
bn_std = np.std(bn_all, axis=0)

# Try to load predictions if they exist
try:
    # You can modify this to load your actual predictions
    from pysr_symbolic_regression import model, X
    y_pred = model.predict(X)
    has_model = True
except:
    has_model = False
    y_pred = None

# Create figure with subplots
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Data overview - all 500 trials
ax1 = axes[0, 0]
for i in range(min(50, bn_all.shape[0])):  # Plot subset for clarity
    ax1.plot(n, bn_all[i], alpha=0.1, color='blue')
ax1.plot(n, bn_avg, 'r-', linewidth=2, label='Average (bn_avg)')
ax1.fill_between(n, bn_avg - bn_std, bn_avg + bn_std, alpha=0.3, color='red', label='+/-1 std')
ax1.set_xlabel('n')
ax1.set_ylabel('bn')
ax1.set_title('bn_all: 500 trials (showing 50) with average')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Averaged data with fit
ax2 = axes[0, 1]
ax2.scatter(n, bn_avg, c='blue', s=30, label='Data (averaged)', zorder=3)
if has_model and y_pred is not None:
    ax2.plot(n, y_pred, 'r-', linewidth=2, label='PySR Fit', zorder=2)
ax2.set_xlabel('n')
ax2.set_ylabel('bn (averaged)')
ax2.set_title('Symbolic Regression Fit')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Residuals
ax3 = axes[1, 0]
if has_model and y_pred is not None:
    residuals = bn_avg - y_pred
    ax3.scatter(n, residuals, c='green', s=30)
    ax3.axhline(y=0, color='r', linestyle='--')
    ax3.set_xlabel('n')
    ax3.set_ylabel('Residuals (data - fit)')
    ax3.set_title('Residuals')
    ax3.grid(True, alpha=0.3)
else:
    ax3.text(0.5, 0.5, 'Run pysr_symbolic_regression.py first', 
             ha='center', va='center', transform=ax3.transAxes)

# Plot 4: Relative error
ax4 = axes[1, 1]
if has_model and y_pred is not None:
    rel_error = np.abs(residuals) / np.abs(bn_avg) * 100
    ax4.scatter(n, rel_error, c='purple', s=30)
    ax4.set_xlabel('n')
    ax4.set_ylabel('Relative Error (%)')
    ax4.set_title('Relative Error')
    ax4.grid(True, alpha=0.3)
else:
    ax4.text(0.5, 0.5, 'Run pysr_symbolic_regression.py first', 
             ha='center', va='center', transform=ax4.transAxes)

plt.tight_layout()
plt.savefig('pysr_visualization.png', dpi=150, bbox_inches='tight')
print("Visualization saved to 'pysr_visualization.png'")
plt.show()
