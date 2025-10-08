import numpy as np

def generate_cost_matrix(n: int, seed: int = None):
    if seed is not None:
        np.random.seed(seed)
    
    points = np.random.rand(n, 2)
    
    cost_matrix = np.sqrt(((points[:, None, :] - points[None, :, :]) ** 2).sum(axis=2))
    
    np.fill_diagonal(cost_matrix, 0)
    
    return cost_matrix

if __name__ == "__main__":
    n = 30
    output_file = "data/cost_matrix.txt"
    
    C = generate_cost_matrix(n)
    
    # Save as space-separated text
    np.savetxt(output_file, C, fmt="%.4f")
    
    print(f"Generated {n}x{n} cost matrix in {output_file}")