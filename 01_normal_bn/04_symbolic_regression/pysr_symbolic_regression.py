"""
Symbolic Regression using PySR to find bn(n) dependence.
Averages bn_all (500 trials) and fits using +, -, *, /, inv(), sqr(), sqrt() operators.
"""

import numpy as np
import scipy.io
from pysr import PySRRegressor

# Load the .mat file
mat_data = scipy.io.loadmat('optimal_bn_500times.mat')

# Extract variables
n = mat_data['n'].flatten()  # Shape: (99,)
bn_all = mat_data['bn_all']  # Shape: (500, 99)

# Average bn_all over 500 trials to get a single vector
bn_avg = np.mean(bn_all, axis=0)  # Shape: (99,)

print(f"n shape: {n.shape}")
print(f"bn_avg shape: {bn_avg.shape}")
print(f"n range: [{n.min()}, {n.max()}]")
print(f"bn_avg range: [{bn_avg.min():.6f}, {bn_avg.max():.6f}]")

# Prepare data for PySR (reshape to column vectors)
X = n.reshape(-1, 1)  # Features: n
y = bn_avg.flatten()  # Target: averaged bn values

# Define the operators for LFF (Linear Fractional Function) form
# LFF form: (a + b*n) / (c + d*n) or similar rational functions
# Using: +, -, *, /, inv (inverse), sqr (square), sqrt (square root)

model = PySRRegressor(
    niterations=1000,
    populations=100,
    population_size=300,
    maxsize=30,                # max expression complexity
    maxdepth=6,                # max expression depth
    
    # Binary operators: +, -, *, /
    binary_operators=[
        "+",
        "-",
        "*",
        "/",
    ],
    
    # Unary operators: inv(x) = 1/x, sqr(x) = x^2, sqrt(x) = safe_sqrt(x)
    unary_operators=[
        "inv(x) = 1/x",
        "sqr(x) = x^2",
        "safe_sqrt(x) = x >= 0 ? sqrt(x) : convert(typeof(x), NaN)",
    ],
    
    # Extra symbolic mappings (optional, for pretty printing)
    extra_sympy_mappings={
        "inv": lambda x: 1/x,
        "sqr": lambda x: x**2,
        "safe_sqrt": lambda x: x**0.5,
    },
    
    # Loss function: mean squared error
    elementwise_loss="loss(x, y) = (x - y)^2",
    
    # Complexity penalties
    parsimony=0.01,
    
    # Progress bar
    progress=True,
    
    # Random seed for reproducibility
    random_state=42,
    
    # Multiprocessing
    procs=4,
    
    # Output
    verbosity=1,
)

print("\n" + "="*60)
print("Starting Symbolic Regression with PySR")
print("="*60)
print(f"Operators: +, -, *, /, inv(), sqr(), safe_sqrt()")
print(f"Target form: LFF (Linear Fractional Function)")
print()

# Fit the model
model.fit(X, y)

print("\n" + "="*60)
print("Results")
print("="*60)

# Display all equations
print("\nAll discovered equations (sorted by score):")
print(model)

# Get the best equation
print("\n" + "-"*60)
print("Best equation (lowest loss):")
print(f"  Equation: {model.get_best()}")
print(f"  SymPy: {model.sympy()}")
print(f"  LaTeX: {model.latex()}")

# Predictions
y_pred = model.predict(X)

# Calculate metrics
mse = np.mean((y - y_pred)**2)
rmse = np.sqrt(mse)
mae = np.mean(np.abs(y - y_pred))

print(f"\nMetrics:")
print(f"  MSE:  {mse:.10f}")
print(f"  RMSE: {rmse:.10f}")
print(f"  MAE:  {mae:.10f}")

# Save results
results = {
    'n': n,
    'bn_avg': bn_avg,
    'bn_pred': y_pred,
    'best_equation': str(model.sympy()),
    'latex_equation': model.latex(),
    'mse': mse,
    'rmse': rmse,
    'mae': mae,
}

# Save to file
with open('pysr_results.txt', 'w') as f:
    f.write("PySR Symbolic Regression Results\n")
    f.write("="*60 + "\n\n")
    f.write(f"Operators: +, -, *, /, inv(), sqr(), safe_sqrt()\n")
    f.write(f"Data points: {len(n)}\n")
    f.write(f"Trials averaged: 500\n\n")
    f.write("Best Equation:\n")
    f.write(f"  SymPy: {model.sympy()}\n")
    f.write(f"  LaTeX: {model.latex()}\n\n")
    f.write(f"Metrics:\n")
    f.write(f"  MSE:  {mse:.10f}\n")
    f.write(f"  RMSE: {rmse:.10f}\n")
    f.write(f"  MAE:  {mae:.10f}\n\n")
    f.write("Top 10 Equations:\n")
    f.write(str(model))

print("\n" + "="*60)
print("Results saved to 'pysr_results.txt'")
print("="*60)
