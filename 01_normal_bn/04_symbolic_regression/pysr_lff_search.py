"""
PySR Symbolic Regression focused on LFF (Linear Fractional Function) form.

LFF form: bn(n) = (a + b*n) / (c + d*n)  or similar rational functions

This script uses constrained search to favor rational function forms.
"""

import numpy as np
import scipy.io
from pysr import PySRRegressor

# Load the .mat file
mat_data = scipy.io.loadmat('optimal_bn_500times.mat')

# Extract and prepare data
n = mat_data['n'].flatten()  # Shape: (99,)
bn_all = mat_data['bn_all']  # Shape: (500, 99)
bn_avg = np.mean(bn_all, axis=0)

# Also check if bn_mean is already computed (should be same as our average)
bn_mean_file = mat_data.get('bn_mean', None)
if bn_mean_file is not None:
    bn_mean_file = bn_mean_file.flatten()
    print(f"bn_mean from file vs computed: max diff = {np.max(np.abs(bn_mean_file - bn_avg)):.2e}")
    # Use the file's bn_mean if it's essentially the same
    bn_avg = bn_mean_file

print(f"Data shape: n={n.shape}, bn_avg={bn_avg.shape}")
print(f"n range: [{n.min()}, {n.max()}]")
print(f"bn range: [{bn_avg.min():.6f}, {bn_avg.max():.6f}]")

# Prepare data
X = n.reshape(-1, 1)
y = bn_avg

# LFF-focused model configuration
# We want to find forms like: (a + b*n) / (c + d*n)
model = PySRRegressor(
    niterations=200,           # More iterations for thorough search
    populations=40,
    population_size=100,
    
    # Limit size to encourage simpler rational forms
    maxsize=15,
    maxdepth=4,
    
    # Binary operators
    binary_operators=["+", "-", "*", "/"],
    
    # Unary operators - inv() is key for LFF forms
    unary_operators=["inv(x) = 1/x"],
    
    # Nested constraints to encourage rational function structure
    # Allow division at top level to get (num)/(denom) form
    nested_constraints={
        "/": {"/": 0},  # Prevent nested divisions like a/(b/c)
    },
    
    # Complexity of each operator
    complexity_of_operators={"/": 2, "+": 1, "-": 1, "*": 1},
    complexity_of_constants=1,
    
    # Parsimony to favor simpler expressions
    parsimony=0.05,
    
    # Loss function (elementwise)
    elementwise_loss="loss(x, y) = (x - y)^2",
    
    # Annealing schedule
    alpha=0.1,
    
    # Seed
    random_state=42,
    
    # Parallel processing
    procs=4,
    
    # Output
    progress=True,
    verbosity=1,
    
    # Save state
    temp_equation_file=True,
    delete_tempfiles=False,
    
    # Extra sympy mappings
    extra_sympy_mappings={"inv": lambda x: 1/x},
    extra_torch_mappings={"inv": lambda x: 1/x},
)

print("\n" + "="*70)
print("PySR Symbolic Regression - LFF Form Search")
print("="*70)
print("Target form: bn(n) = (a + b*n) / (c + d*n)  or similar")
print("Operators: +, -, *, /, inv()")
print()

# Fit
model.fit(X, y)

# Results
print("\n" + "="*70)
print("RESULTS")
print("="*70)

# Get all equations
print("\nAll discovered equations:")
print(model)

# Best equation
best_eq = model.get_best()
print(f"\nBest equation:")
print(f"  PySR: {best_eq}")
print(f"  SymPy: {model.sympy()}")
print(f"  LaTeX: ${model.latex()}$")

# Predictions and metrics
y_pred = model.predict(X)
mse = np.mean((y - y_pred)**2)
rmse = np.sqrt(mse)
mape = np.mean(np.abs((y - y_pred) / y)) * 100
r2 = 1 - np.sum((y - y_pred)**2) / np.sum((y - np.mean(y))**2)

print(f"\nMetrics:")
print(f"  MSE:   {mse:.10e}")
print(f"  RMSE:  {rmse:.10e}")
print(f"  MAPE:  {mape:.6f}%")
print(f"  R^2:    {r2:.10f}")

# Save results
with open('lff_results.txt', 'w') as f:
    f.write("="*70 + "\n")
    f.write("PySR LFF Symbolic Regression Results\n")
    f.write("="*70 + "\n\n")
    
    f.write("Configuration:\n")
    f.write(f"  Iterations: 200\n")
    f.write(f"  Populations: 40\n")
    f.write(f"  Operators: +, -, *, /, inv()\n")
    f.write(f"  Max size: 15\n")
    f.write(f"  Max depth: 4\n\n")
    
    f.write("Data:\n")
    f.write(f"  Points: {len(n)}\n")
    f.write(f"  Trials averaged: 500\n")
    f.write(f"  n range: [{n.min()}, {n.max()}]\n\n")
    
    f.write("Best Equation:\n")
    f.write(f"  SymPy: {model.sympy()}\n")
    f.write(f"  LaTeX: {model.latex()}\n\n")
    
    f.write("Metrics:\n")
    f.write(f"  MSE:   {mse:.10e}\n")
    f.write(f"  RMSE:  {rmse:.10e}\n")
    f.write(f"  MAPE:  {mape:.6f}%\n")
    f.write(f"  R^2:    {r2:.10f}\n\n")
    
    f.write("Top Equations (by score):\n")
    f.write(str(model))

print("\n" + "="*70)
print("Results saved to 'lff_results.txt' and 'lff_equations.csv'")
print("="*70)

# Optional: Quick visualization
print("\nGenerating quick visualization...")
try:
    import matplotlib.pyplot as plt
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Fit plot
    ax1 = axes[0]
    ax1.scatter(n, y, c='blue', s=40, label='Data (avg of 500)', zorder=3, alpha=0.7)
    ax1.plot(n, y_pred, 'r-', linewidth=2, label='PySR Fit', zorder=2)
    ax1.set_xlabel('n')
    ax1.set_ylabel('bn')
    ax1.set_title('Symbolic Regression Fit (LFF Search)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Residuals
    ax2 = axes[1]
    residuals = y - y_pred
    ax2.scatter(n, residuals, c='green', s=40)
    ax2.axhline(y=0, color='r', linestyle='--')
    ax2.set_xlabel('n')
    ax2.set_ylabel('Residuals')
    ax2.set_title('Residuals (Data - Fit)')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('lff_fit.png', dpi=150, bbox_inches='tight')
    print("Visualization saved to 'lff_fit.png'")
except Exception as e:
    print(f"Could not generate visualization: {e}")
