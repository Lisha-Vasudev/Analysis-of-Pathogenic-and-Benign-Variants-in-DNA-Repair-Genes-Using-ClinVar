#!/usr/bin/env python
"""
Orchestration script for variant analysis workflow.

This script documents and can optionally run the complete analysis pipeline:
1. Data preparation (Notebook 01)
2. Analysis and interpretation (Notebook 02)

For interactive exploration, run notebooks manually with:
    jupyter notebook notebooks/01_data_preparation.ipynb
    jupyter notebook notebooks/02_analysis_and_interpretation.ipynb

For batch execution, ensure all dependencies are installed and run:
    python run_analysis.py
"""

import os
import sys
import subprocess
from pathlib import Path


def get_project_root():
    """Get the root directory of the project."""
    return Path(__file__).parent.resolve()


def check_dependencies():
    """Check if all required packages are installed."""
    required_packages = [
        'pandas',
        'numpy',
        'matplotlib',
        'seaborn',
        'jupyter',
        'ipython',
    ]
    
    missing = []
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing.append(package)
    
    if missing:
        print(f"ERROR: Missing required packages: {', '.join(missing)}")
        print(f"Install with: pip install -r requirements.txt")
        return False
    return True


def check_data_setup():
    """Check if data directories are properly configured."""
    project_root = get_project_root()
    
    required_dirs = [
        'data/raw',
        'data/processed',
        'notebooks',
        'scripts',
        'results/figures',
        'results/tables',
        'docs'
    ]
    
    for dir_path in required_dirs:
        full_path = project_root / dir_path
        if not full_path.exists():
            print(f"Creating directory: {dir_path}")
            full_path.mkdir(parents=True, exist_ok=True)
    
    return True


def run_notebooks_sequentially():
    """Run both notebooks in sequence using jupyter nbconvert."""
    project_root = get_project_root()
    
    notebooks = [
        'notebooks/01_data_preparation.ipynb',
        'notebooks/02_analysis_and_interpretation.ipynb'
    ]
    
    print("\n" + "="*70)
    print("VARIANT ANALYSIS PIPELINE")
    print("="*70)
    print(f"Project root: {project_root}")
    
    for i, notebook in enumerate(notebooks, 1):
        notebook_path = project_root / notebook
        
        print(f"\n[{i}/{len(notebooks)}] Running: {notebook}")
        print("-" * 70)
        
        if not notebook_path.exists():
            print(f"ERROR: Notebook not found at {notebook_path}")
            return False
        
        try:
            # Use nbconvert to run the notebook
            output_path = notebook_path.with_name(
                notebook_path.stem + '_executed.ipynb'
            )
            
            cmd = [
                'jupyter', 'nbconvert',
                '--to', 'notebook',
                '--execute',
                '--ExecutePreprocessor.timeout=600',
                '--ExecutePreprocessor.kernel_name=python3',
                f'--output={output_path}',
                str(notebook_path)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"ERROR: Notebook execution failed")
                print(f"stdout: {result.stdout}")
                print(f"stderr: {result.stderr}")
                return False
            
            print(f"✓ Completed: {notebook}")
            print(f"  Output saved to: {output_path}")
            
        except FileNotFoundError:
            print("ERROR: jupyter nbconvert not found.")
            print("Install with: pip install jupyterlab nbconvert")
            return False
    
    return True


def print_workflow_instructions():
    """Print instructions for running the analysis."""
    project_root = get_project_root()
    
    instructions = f"""
{'='*70}
VARIANT ANALYSIS WORKFLOW - SETUP COMPLETE
{'='*70}

This project analyzes pathogenic and benign variants in DNA repair genes
from the ClinVar public database.

PROJECT STRUCTURE:
  notebooks/          - Jupyter notebooks for analysis
  data/               - Raw and processed data
  scripts/            - Reusable Python modules
  results/            - Output figures and tables
  docs/               - Methods and documentation

QUICK START (RECOMMENDED - Interactive):
  1. Open a terminal and navigate to: {project_root}
  2. Start Jupyter:
     jupyter notebook
  3. Open and run: notebooks/01_data_preparation.ipynb
  4. Then run: notebooks/02_analysis_and_interpretation.ipynb

AUTOMATED EXECUTION (Batch mode):
  python run_analysis.py --execute
  (Warning: downloads ~700 MB ClinVar data on first run)

DEPENDENCIES:
  All required packages are listed in requirements.txt
  Install with: pip install -r requirements.txt

DATA:
  ClinVar data is automatically downloaded on first run (Notebook 01).
  Raw data saved to: data/raw/variant_summary.txt
  Processed data saved to: data/processed/variants_filtered.csv

OUTPUTS:
  Figures: results/figures/ (PNG, 300 DPI)
  Tables:  results/tables/  (CSV format)

DOCUMENTATION:
  - README.md         - Project overview
  - docs/methods.md   - Detailed methods and interpretation
  - Notebook cells    - Inline documentation and explanations

FOR QUESTIONS:
  See docs/methods.md for detailed methods, limitations, and interpretation.

{'='*70}
"""
    
    print(instructions)


def main():
    """Main orchestration function."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Variant Analysis Workflow Orchestration'
    )
    parser.add_argument(
        '--execute',
        action='store_true',
        help='Execute notebooks in sequence (batch mode)'
    )
    parser.add_argument(
        '--setup',
        action='store_true',
        help='Setup directories and check dependencies'
    )
    
    args = parser.parse_args()
    
    # Always check setup
    print("Checking project setup...")
    if not check_dependencies():
        sys.exit(1)
    if not check_data_setup():
        sys.exit(1)
    print("✓ Setup check passed\n")
    
    # Execute if requested
    if args.execute:
        print("Starting notebook execution...")
        if not run_notebooks_sequentially():
            sys.exit(1)
        print("\n✓ Analysis complete!")
        print("Results saved to results/ directory")
    else:
        # Print workflow instructions
        print_workflow_instructions()


if __name__ == '__main__':
    main()
    print("  Results:")
    print("    - results/variant_distributions.png")
    print("    - results/chromosome_distribution.png")
    print("    - results/clinical_significance.png")
    print("    - results/position_density.png")
    print("\nNext steps:")
    print("  1. Open notebooks/analysis.ipynb for detailed analysis")
    print("  2. Modify load_data.py to use your own dataset")
    print("  3. Customize filtering thresholds in clean_data.py")


if __name__ == '__main__':
    main()
