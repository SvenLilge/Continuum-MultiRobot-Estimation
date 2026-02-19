#!/usr/bin/env python3
"""
Script to plot estimated disk poses against ground truth data.
Generates comparison plots for position and rotation of all 7 disks.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from scipy.linalg import logm


def load_data(csv_file):
    """Load CSV file and return dataframe."""
    df = pd.read_csv(csv_file)
    return df


def euler_to_rotation_matrix(rotX, rotY, rotZ):
    """Convert Euler angles to rotation matrix using XYZ convention.
    
    Args:
        rotX, rotY, rotZ: Euler angles in radians (XYZ convention)
    
    Returns:
        3x3 rotation matrix R = Rx * Ry * Rz
    """
    cx, sx = np.cos(rotX), np.sin(rotX)
    cy, sy = np.cos(rotY), np.sin(rotY)
    cz, sz = np.cos(rotZ), np.sin(rotZ)
    
    Rx = np.array([[1,  0,   0],
                   [0, cx, -sx],
                   [0, sx,  cx]])
    
    Ry = np.array([[cy, 0, sy],
                   [0,  1,  0],
                   [-sy, 0, cy]])
    
    Rz = np.array([[cz, -sz, 0],
                   [sz,  cz, 0],
                   [0,   0,  1]])
    
    return Rx @ Ry @ Rz


def compute_orientation_error(R1, R2):
    """Compute orientation error between two rotation matrices.
    
    Args:
        R1, R2: 3x3 rotation matrices
    
    Returns:
        Scalar error: norm of log(inv(R1) * R2)
    """
    R_error = np.linalg.inv(R1) @ R2
    log_R_error = logm(R_error)
    # Frobenius norm of the matrix logarithm
    return np.linalg.norm(log_R_error)


def plot_disk_comparison(ax, gt_data, est_data, disk_idx, component, ylabel):
    """Plot a single disk component (position or rotation) comparison."""
    gt_col = f'disk_{disk_idx}_{component}'
    est_col = f'disk_{disk_idx}_{component}'
    
    t_gt = gt_data['timestamp'].values
    t_est = est_data['timestamp'].values
    
    gt_vals = gt_data[gt_col].values if gt_col in gt_data.columns else None
    est_vals = est_data[est_col].values if est_col in est_data.columns else None
    
    if gt_vals is not None:
        ax.plot(t_gt, gt_vals, 'b-', linewidth=2, label='Ground Truth', alpha=0.7)
    if est_vals is not None:
        ax.plot(t_est, est_vals, 'r--', linewidth=2, label='Estimated', alpha=0.7)
    
    ax.set_ylabel(ylabel, fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=8)


def plot_all_disks(gt_data, est_data, output_dir=None):
    """Generate plots for all 7 disks with position and rotation data."""
    
    n_disks = 7
    components = ['x', 'y', 'z']  # position components
    rot_components = ['rotX', 'rotY', 'rotZ']  # rotation components
    
    # Create figure with subplots for each disk (6 subplots per disk: 3 pos + 3 rot)
    fig, axes = plt.subplots(n_disks, 6, figsize=(18, 14))
    fig.suptitle('Disk Pose Estimation: Ground Truth vs Estimated', fontsize=16, fontweight='bold')
    
    # Plot position and rotation for each disk
    for disk_idx in range(n_disks):
        # Position plots
        for comp_idx, comp in enumerate(components):
            ax = axes[disk_idx, comp_idx]
            plot_disk_comparison(ax, gt_data, est_data, disk_idx, comp, f'Disk {disk_idx}\n{comp} (m)')
            if disk_idx == 0 and comp_idx == 0:
                ax.legend(loc='upper right', fontsize=8)
            else:
                ax.legend([])
        
        # Rotation plots
        for comp_idx, comp in enumerate(rot_components):
            ax = axes[disk_idx, 3 + comp_idx]
            plot_disk_comparison(ax, gt_data, est_data, disk_idx, comp, f'Disk {disk_idx}\n{comp} (rad)')
            ax.legend([])
        
        # Add disk label on the left
        fig.text(0.02, 0.95 - disk_idx * (1.0 / n_disks), f'Disk {disk_idx}', 
                fontsize=10, fontweight='bold', va='center')
    
    fig.text(0.5, 0.02, 'Time (s)', ha='center', fontsize=12, fontweight='bold')
    
    plt.tight_layout(rect=[0.03, 0.03, 1, 0.97])
    
    # Save figure if output directory is specified
    if output_dir:
        output_path = Path(output_dir) / 'disk_comparison.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved to {output_path}")
    
    plt.show()


def plot_individual_disks(gt_data, est_data, output_dir=None):
    """Generate individual plots for each disk with position, orientation, and error subplots."""
    
    n_disks = 7
    components = ['x', 'y', 'z']
    rot_components = ['rotX', 'rotY', 'rotZ']
    
    for disk_idx in range(n_disks):
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        fig.suptitle(f'Disk {disk_idx} - Pose Comparison (Ground Truth vs Estimated)', 
                    fontsize=14, fontweight='bold')
        
        t_gt = gt_data['timestamp'].values
        t_est = est_data['timestamp'].values
        
        # Position subplot (top-left)
        ax_pos = axes[0, 0]
        for comp_idx, comp in enumerate(components):
            # GT column format: "disk_0 x (m)"
            gt_col = f'disk_{disk_idx} {comp} (m)'
            # EST column format: "disk_0_x"
            est_col = f'disk_{disk_idx}_{comp}'
            
            if gt_col in gt_data.columns:
                ax_pos.plot(t_gt, gt_data[gt_col].values, 'o-', linewidth=2, 
                           markersize=4, label=f'GT {comp}', alpha=0.7)
            
            if est_col in est_data.columns:
                ax_pos.plot(t_est, est_data[est_col].values, 's--', linewidth=2, 
                           markersize=3, label=f'Est {comp}', alpha=0.7)
        
        ax_pos.set_xlabel('Time (s)', fontsize=11)
        ax_pos.set_ylabel('Position (m)', fontsize=11)
        ax_pos.set_title('Position (X, Y, Z)', fontsize=12, fontweight='bold')
        ax_pos.grid(True, alpha=0.3)
        ax_pos.legend(loc='best', fontsize=9, ncol=2)
        
        # Orientation subplot (top-right)
        ax_rot = axes[0, 1]
        for comp_idx, comp in enumerate(rot_components):
            # GT column format: "disk_0 rotX (rad)"
            gt_col = f'disk_{disk_idx} {comp} (rad)'
            # EST column format: "disk_0_rotX"
            est_col = f'disk_{disk_idx}_{comp}'
            
            if gt_col in gt_data.columns:
                ax_rot.plot(t_gt, gt_data[gt_col].values, 'o-', linewidth=2, 
                           markersize=4, label=f'GT {comp}', alpha=0.7)
            
            if est_col in est_data.columns:
                ax_rot.plot(t_est, est_data[est_col].values, 's--', linewidth=2, 
                           markersize=3, label=f'Est {comp}', alpha=0.7)
        
        ax_rot.set_xlabel('Time (s)', fontsize=11)
        ax_rot.set_ylabel('Rotation (rad)', fontsize=11)
        ax_rot.set_title('Orientation (rotX, rotY, rotZ)', fontsize=12, fontweight='bold')
        ax_rot.grid(True, alpha=0.3)
        ax_rot.legend(loc='best', fontsize=9, ncol=2)
        
        # Position error subplot (bottom-left)
        ax_pos_err = axes[1, 0]
        pos_errors = []
        pos_error_times = []
        
        for i, t in enumerate(t_est):
            # Find closest GT timestamp
            idx_gt = np.argmin(np.abs(t_gt - t))
            
            # Get GT position
            gt_x_col = f'disk_{disk_idx} x (m)'
            gt_y_col = f'disk_{disk_idx} y (m)'
            gt_z_col = f'disk_{disk_idx} z (m)'
            
            # Get estimated position
            est_x_col = f'disk_{disk_idx}_x'
            est_y_col = f'disk_{disk_idx}_y'
            est_z_col = f'disk_{disk_idx}_z'
            
            if all(col in gt_data.columns for col in [gt_x_col, gt_y_col, gt_z_col]) and \
               all(col in est_data.columns for col in [est_x_col, est_y_col, est_z_col]):
                
                gt_pos = np.array([
                    gt_data[gt_x_col].iloc[idx_gt],
                    gt_data[gt_y_col].iloc[idx_gt],
                    gt_data[gt_z_col].iloc[idx_gt]
                ])
                
                est_pos = np.array([
                    est_data[est_x_col].iloc[i],
                    est_data[est_y_col].iloc[i],
                    est_data[est_z_col].iloc[i]
                ])
                
                pos_error = np.linalg.norm(gt_pos - est_pos)
                pos_errors.append(pos_error)
                pos_error_times.append(t)
        
        if pos_errors:
            ax_pos_err.plot(pos_error_times, pos_errors, 'r-', linewidth=2, label='Position Error')
            ax_pos_err.set_xlabel('Time (s)', fontsize=11)
            ax_pos_err.set_ylabel('Error (m)', fontsize=11)
            ax_pos_err.set_title('Position Error Norm', fontsize=12, fontweight='bold')
            ax_pos_err.grid(True, alpha=0.3)
            ax_pos_err.legend(loc='best', fontsize=9)
        
        # Orientation error subplot (bottom-right)
        ax_rot_err = axes[1, 1]
        rot_errors = []
        rot_error_times = []
        
        for i, t in enumerate(t_est):
            # Find closest GT timestamp
            idx_gt = np.argmin(np.abs(t_gt - t))
            
            # Get GT Euler angles
            gt_rotX_col = f'disk_{disk_idx} rotX (rad)'
            gt_rotY_col = f'disk_{disk_idx} rotY (rad)'
            gt_rotZ_col = f'disk_{disk_idx} rotZ (rad)'
            
            # Get estimated Euler angles
            est_rotX_col = f'disk_{disk_idx}_rotX'
            est_rotY_col = f'disk_{disk_idx}_rotY'
            est_rotZ_col = f'disk_{disk_idx}_rotZ'
            
            if all(col in gt_data.columns for col in [gt_rotX_col, gt_rotY_col, gt_rotZ_col]) and \
               all(col in est_data.columns for col in [est_rotX_col, est_rotY_col, est_rotZ_col]):
                
                # Construct GT rotation matrix
                gt_rotX = gt_data[gt_rotX_col].iloc[idx_gt]
                gt_rotY = gt_data[gt_rotY_col].iloc[idx_gt]
                gt_rotZ = gt_data[gt_rotZ_col].iloc[idx_gt]
                R_gt = euler_to_rotation_matrix(gt_rotX, gt_rotY, gt_rotZ)
                
                # Construct estimated rotation matrix
                est_rotX = est_data[est_rotX_col].iloc[i]
                est_rotY = est_data[est_rotY_col].iloc[i]
                est_rotZ = est_data[est_rotZ_col].iloc[i]
                R_est = euler_to_rotation_matrix(est_rotX, est_rotY, est_rotZ)
                
                # Compute orientation error
                rot_error = compute_orientation_error(R_est, R_gt)
                rot_errors.append(rot_error)
                rot_error_times.append(t)
        
        if rot_errors:
            ax_rot_err.plot(rot_error_times, rot_errors, 'b-', linewidth=2, label='Orientation Error')
            ax_rot_err.set_xlabel('Time (s)', fontsize=11)
            ax_rot_err.set_ylabel('Error (rad)', fontsize=11)
            ax_rot_err.set_title('Orientation Error (||log(R)||)', fontsize=12, fontweight='bold')
            ax_rot_err.grid(True, alpha=0.3)
            ax_rot_err.legend(loc='best', fontsize=9)
        
        plt.tight_layout()
        
        # Save figure if output directory is specified
        if output_dir:
            output_path = Path(output_dir) / f'disk_{disk_idx}_comparison.png'
            plt.savefig(output_path, dpi=150, bbox_inches='tight')
            print(f"Figure saved to {output_path}")
        
        plt.show()


def compute_statistics(gt_data, est_data):
    """Compute and print statistics for the comparison."""
    
    n_disks = 7
    components = ['x', 'y', 'z']
    rot_components = ['rotX', 'rotY', 'rotZ']
    
    print("\n" + "="*80)
    print("ESTIMATION STATISTICS")
    print("="*80)
    
    for disk_idx in range(n_disks):
        print(f"\nDisk {disk_idx}:")
        print("-" * 40)
        
        # Position statistics
        print("  Position (meters):")
        for comp in components:
            gt_col = f'disk_{disk_idx}_{comp}'
            if gt_col in gt_data.columns and gt_col in est_data.columns:
                # Interpolate est_data to match gt_data timestamps
                est_interp = np.interp(gt_data['timestamp'], est_data['timestamp'], 
                                      est_data[gt_col])
                error = gt_data[gt_col].values - est_interp
                mae = np.mean(np.abs(error))
                rmse = np.sqrt(np.mean(error**2))
                print(f"    {comp}: MAE={mae:.6f}m, RMSE={rmse:.6f}m")
        
        # Rotation statistics
        print("  Rotation (radians):")
        for comp in rot_components:
            gt_col = f'disk_{disk_idx}_{comp}'
            if gt_col in gt_data.columns and gt_col in est_data.columns:
                # Interpolate est_data to match gt_data timestamps
                est_interp = np.interp(gt_data['timestamp'], est_data['timestamp'], 
                                      est_data[gt_col])
                error = gt_data[gt_col].values - est_interp
                mae = np.mean(np.abs(error))
                rmse = np.sqrt(np.mean(error**2))
                print(f"    {comp}: MAE={mae:.6f}rad, RMSE={rmse:.6f}rad")


def main():
    parser = argparse.ArgumentParser(
        description='Plot estimated disk poses against ground truth data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python plot_estimates.py
  python plot_estimates.py --individual
  python plot_estimates.py --stats
  python plot_estimates.py --gt custom_gt.csv --est custom_est.csv
        """
    )
    
    parser.add_argument('--gt', type=str, default='../data/RSS2026/base 1/multiCR oscillation 1/dataVicon.csv',
                       help='Path to ground truth CSV file (default: ../data/RSS2026/base 1/multiCR oscillation 1/dataVicon.csv)')
    parser.add_argument('--est', type=str, default='../data/RSS2026/base 1/multiCR oscillation 1/GP_estimates.csv',
                       help='Path to estimated poses CSV file (default: ../data/RSS2026/base 1/multiCR oscillation 1/GP_estimates.csv)')
    parser.add_argument('--individual', action='store_true',
                       help='Plot individual figures for each disk instead of combined plot')
    parser.add_argument('--stats', action='store_true',
                       help='Compute and print statistics')
    parser.add_argument('--output', type=str, default=None,
                       help='Output directory for saving plots (default: no saving)')
    
    args = parser.parse_args()
    
    # Check if files exist
    gt_path = Path(args.gt)
    est_path = Path(args.est)
    
    if not gt_path.exists():
        print(f"Error: Ground truth file not found: {gt_path}")
        return 1
    
    if not est_path.exists():
        print(f"Error: Estimated poses file not found: {est_path}")
        return 1
    
    print("Loading data...")
    gt_data = load_data(gt_path)
    est_data = load_data(est_path)
    
    print(f"Ground truth data shape: {gt_data.shape}")
    print(f"Estimated data shape: {est_data.shape}")
    
    # Create output directory if specified
    if args.output:
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output directory: {output_dir}")
    else:
        output_dir = None
    
    # Generate plots
    if args.individual:
        print("Generating individual disk plots...")
        plot_individual_disks(gt_data, est_data, output_dir)
    else:
        print("Generating combined disk plot...")
        plot_all_disks(gt_data, est_data, output_dir)
    
    # Compute statistics if requested
    if args.stats:
        compute_statistics(gt_data, est_data)
    
    print("Done!")
    return 0


if __name__ == '__main__':
    exit(main())
