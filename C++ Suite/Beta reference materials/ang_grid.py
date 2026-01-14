import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

def repulsion_orientations(n_orientations, gamma_fixed=0.0, C=1e-3, tol=1e-4, max_iter=5000, seed=None):
    def fibonacci_sphere(N, perturb_strength=0.1):
        golden_ratio = (1 + 5**0.5) / 2
        i = np.arange(N)
        z = 1 - 2*(i + 0.5)/N
        theta = np.arccos(z)
        phi = 2 * np.pi * i / golden_ratio
        phi = np.mod(phi, 2*np.pi)
    
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = z

        points = np.stack((x, y, z), axis=1)

        # Add random perturbation perpendicular to each point
        perturbed = []
        for p in points:
            # Generate random vector
            rand = np.random.randn(3)
            # Project onto tangent plane
            rand -= np.dot(rand, p) * p
            rand = rand / np.linalg.norm(rand)
            # Add small perturbation
            p_new = p + perturb_strength * rand
            p_new /= np.linalg.norm(p_new)
            perturbed.append(p_new)
    
        return np.array(perturbed)

    OR = fibonacci_sphere(n_orientations)

    for iteration in range(max_iter):
        displacements = np.zeros_like(OR)
        max_disp = 0

        for i in range(n_orientations):
            for j in range(i + 1, n_orientations):
                vi = OR[i]
                vj = OR[j]
                dot = np.clip(np.dot(vi, vj), -1.0, 1.0)
                theta_ij = np.arccos(dot)

                V = np.cross(vj, vi)
                dP_i = np.cross(V, vi)
                dP_j = np.cross(vj, V)

                dP_i /= np.linalg.norm(dP_i)
                dP_j /= np.linalg.norm(dP_j)

                F = C / (theta_ij)**2

                displacements[i] += F * dP_i
                displacements[j] += F * dP_j

        for i in range(n_orientations):
            OR[i] += displacements[i]
            OR[i] /= np.linalg.norm(OR[i])
            max_disp = max(max_disp, np.linalg.norm(displacements[i]))

        if max_disp < tol:
            print(f"Converged after {iteration + 1} iterations.")
            break
    else:
        print(f"Did not converge after {max_iter} iterations.")

    # Convert to Euler angles (ZYZ convention: α, β, γ), with γ fixed
    alpha = np.arctan2(OR[:, 1], OR[:, 0]) + np.pi
    beta = np.arccos(OR[:, 2])
    gamma = np.full_like(alpha, gamma_fixed)

    angle_grid = np.stack((alpha, beta, gamma), axis=1)
    return angle_grid

def visualize_orientations(orientations):
    alpha = np.degrees(orientations[:, 0]) 
    beta = np.degrees(orientations[:, 1])

    a_rad = np.radians(alpha)
    b_rad = np.radians(beta)
    x = np.sin(b_rad) * np.cos(a_rad)
    y = np.sin(b_rad) * np.sin(a_rad)
    z = np.cos(b_rad)
    points = np.vstack((x, y, z)).T

    hull = ConvexHull(points)

    fig = plt.figure(figsize=(12, 6))

    ax1 = fig.add_subplot(1, 2, 1)
    ax1.scatter(alpha, beta, color='black', s=8)
    ax1.set_xlim(0, 360)
    ax1.set_ylim(0, 180)
    ax1.set_xlabel(r'$\alpha$ (deg)')
    ax1.set_ylabel(r'$\beta$ (deg)')
    ax1.set_title('Euler Angles (α, β)')

    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    for simplex in hull.simplices:
        ax2.plot(points[simplex, 0], points[simplex, 1], points[simplex, 2], color='black', linewidth=0.5)

    ax2.set_box_aspect([1, 1, 1])
    ax2.view_init(elev=20, azim=45)
    ax2.set_title('Triangulated Sphere Surface')
    ax2.axis('off')

    plt.tight_layout()
    plt.show()

# Example usage
if __name__ == "__main__":
    orientations = repulsion_orientations(n_orientations=50, gamma_fixed=0.0)
    
    # Uncomment below to visualize
    visualize_orientations(orientations)
