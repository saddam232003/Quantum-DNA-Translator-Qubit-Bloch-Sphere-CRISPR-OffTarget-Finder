# Quantum DNA Translator – 11‑Qubit Simulator

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey)]()

**A state‑of‑the‑art quantum bioinformatics desktop application**  
Translates DNA sequences into proteins using an 11‑qubit oracle, visualises the underlying quantum state in real‑time, and includes an embedded CRISPR off‑target finder – all in a professional dark‑theme GUI.

Designed and developed by **Dr. Muhammad Saddam Khokhar** – [Quantummind.AI](https://quantummind.ai)

---

## Table of Contents
- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage Guide](#usage-guide)
- [Graphical User Interface](#graphical-user-interface)
- [How It Works (Quantum Background)](#how-it-works-quantum-background)
- [Dependencies](#dependencies)
- [Customisation](#customisation)
- [Performance Notes](#performance-notes)
- [Contributing](#contributing)
- [License](#license)
- [Author & Contact](#author--contact)

---

## Overview

The **Quantum DNA Translator** bridges quantum information science and molecular biology.  
It uses a full **11‑qubit state‑vector simulation** to execute the standard (or alternative) genetic code as a single unitary operation.  
The graphical interface shows **six live Bloch spheres** (one per qubit), a real‑time probability chart, a circuit diagram, and a growing protein sequence.  
You can also type a CRISPR guide RNA and instantly find off‑target sites in your DNA sequence.

The entire application is written in **pure Python** using only `numpy`, `matplotlib`, and `scipy` – no external quantum frameworks required.  
It runs comfortably on a modern laptop (e.g., Lenovo Legion Core i9, 16 GB RAM).

---

## Key Features

- ✅ **11‑qubit exact simulation** – state‑vector of 2048 complex amplitudes
- ✅ **16 built‑in genetic codes** (standard, vertebrate mitochondrial, yeast, etc.)
- ✅ **Real‑time Bloch spheres** – visualise entanglement and superposition
- ✅ **Quantum measurement** – output qubits collapse to an amino acid index
- ✅ **CRISPR off‑target finder** – type a guide RNA directly in the GUI
- ✅ **Display options** – verbose three‑letter codes, spaces, show DNA, strand selection
- ✅ **State viewer pop‑up** – codon‑state amplitudes, reduced density matrix, Bloch vectors
- ✅ **Genetic code heatmap**, **DNA image generator**, **3D protein placeholder**
- ✅ **Report generation** – saves a detailed text report with pop‑up display
- ✅ **Noise model** – optional depolarising noise (p=0.01)
- ✅ **Professional dark‑theme GUI** – no overlapping widgets, perfectly spaced layout

---

## Installation

### Prerequisites
- Python 3.8 or higher
- `pip` package manager

### Clone the repository
```bash
git clone https://github.com/saddam232003/Variational-Quantum-Eigensolver-VQE-for-the-MaxCut-problem.git
cd Variational-Quantum-Eigensolver-VQE-for-the-MaxCut-problem
```

*(The translator script is located in the same repository – look for `quantum_dna_translator.py` or the main Python file.)*

### Install dependencies
```bash
pip install numpy matplotlib scipy
```

That’s it – no additional libraries are needed.

---

## Quick Start

1. Run the script:
   ```bash
   python quantum_dna_translator.py
   ```
2. In the DNA input box, type or paste a DNA sequence (letters `A`, `T`, `G`, `C`), or use the default pre‑filled sequence.
3. Click **Translate** to see the protein appear, while the Bloch spheres and circuit diagram animate.
4. Type a 20‑nt CRISPR guide RNA (e.g., `TAGCTAGCTAGCATGCATGC`) and click **Find Off‑Targets** to search for binding sites.
5. Explore the other tools: **State**, **Heatmap**, **DNA Img**, **3D Prot**, **Display Opts**, **Report**.

---

## Usage Guide

### DNA Translation
- **Select genetic code** from the 4×4 grid (default: Standard).
- Check/uncheck **Verbose**, **Spaces**, **Show DNA** and choose **strand** (`forward`, `reverse`, `both`).
- Click **Translate** (or press Enter in the DNA text box).  
  The translation runs codon‑by‑codon with a short animation.

### CRISPR Off‑Target Search
- Enter a 20‑nucleotide guide RNA in the permanent text box (below the DNA input).
- Click **Find Off‑Targets**.  
  A pop‑up window shows all binding sites with up to 3 mismatches, along with the number of mismatches and positions.

### Auxiliary Tools
| Button | Function |
|--------|----------|
| **State** | Opens a viewer with codon amplitudes, reduced density matrix, Bloch vectors, and output probabilities. |
| **Heatmap** | Displays a colour‑coded codon table for the first two bases. |
| **DNA Img** | Generates a visual representation of the DNA sequence (A=red, T=blue, G=green, C=yellow). |
| **3D Prot** | Random coil 3D scatter of the translated amino acids, coloured by property. |
| **Display Opts** | Pop‑up to change display format (verbosity, spacing, show DNA, strand). |
| **Report** | Creates a text report (`quantum_translation_report.txt`) and shows it in a pop‑up. |
| **Noise checkbox** | Enables depolarising noise during translation. |

---

## Graphical User Interface

![GUI Screenshot](screenshot.png) *(Replace with an actual screenshot)*

The GUI is built with `matplotlib` and uses a strict `GridSpec` layout to prevent overlapping.  
The dark theme reduces eye strain and highlights the cyan, gold, and orange accents.

Layout:
- **Left top** – translation circuit diagram
- **Left bottom** – output qubits measurement probability bar chart
- **Right top** – six 3D Bloch spheres (Qubit 0‑5)
- **Right middle** – DNA input, CRISPR guide input, translate button, protein display
- **Right middle‑bottom** – genetic code grid (4×4), control buttons (State, Heatmap, DNA Img, 3D Prot), display options, report, noise

---

## How It Works (Quantum Background)

### Encoding
Each nucleotide is represented by **2 qubits** (A=00, T=01, G=10, C=11).  
A codon (3 nucleotides) → **6 qubits**.  
The output amino acid is encoded in **5 qubits** (0‑19 + stop).  
Total: **11 qubits**.

### Oracle
A single **permutation matrix** \( U \) of size \( 2^{11} \times 2^{11} \) is built for the chosen genetic code.  
It maps \( |\text{codon}\rangle|0\rangle \rightarrow |\text{codon}\rangle|\text{amino acid}\rangle \).

### Simulation
The state vector evolves exactly under the oracle. After each codon, the **reduced density matrix** of each qubit is computed (partial trace) to obtain the **Bloch vector** for visualisation.

### Measurement
The 5 output qubits are measured to obtain the amino acid index.

---

## Dependencies

| Package     | Version (tested) | Purpose          |
|-------------|------------------|------------------|
| `numpy`     | ≥1.21            | State‑vector arithmetic |
| `matplotlib`| ≥3.5             | GUI, 3D plotting, widgets |
| `scipy`     | ≥1.7             | COBYLA optimiser (used only in VQE version) |

*(Note: `scipy` is imported but not actively used in the current translator; it can be removed without effect.)*

---

## Customisation

- **Number of qubits**: Change `n_qubits` in `QuantumState` (be cautious – memory grows as \(2^n\)).
- **Genetic codes**: Add new dictionaries to `GENETIC_CODES`.
- **Translation speed**: Modify the `plt.pause(0.3)` values inside `process_codon`.
- **Noise level**: Adjust `self.noise_prob` (default 0.01).
- **GUI colours**: Edit the `self.bg_main`, `self.cyan`, etc. variables in `QuantumDNAGui.__init__`.

---

## Performance Notes

- **Memory**: The state vector occupies ~32 KB. The entire application uses < 200 MB.
- **CPU**: Translation of a 72‑base DNA sequence completes in under 5 seconds (including animation pauses).
- **Rendering**: Matplotlib’s interactive mode (`plt.ion()`) keeps the GUI responsive; pop‑ups use `plt.show(block=False)`.

---

## Contributing

Contributions are welcome! Possible improvements:
- Add a real quantum backend (IBM Q, Rigetti) via Qiskit/Cirq.
- Implement a more sophisticated protein 3D model (PDB output).
- Extend the CRISPR search to use quantum pattern matching (Grover’s algorithm).
- Translate the UI to other languages.

Please open an issue or pull request on GitHub.

---

## License

This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.  
You are free to use, modify, and distribute the code for any purpose, provided you retain the original copyright notice.

---

## Author & Contact

**Dr. Muhammad Saddam Khokhar**  
Founder & Lead Researcher – [Quantummind.AI](https://quantummind.ai)  

📧 saddam_khokhar@hotmail.com  

If you use this software in your research or teaching, please cite the author and repository.  
For commercial inquiries or collaboration, feel free to reach out.

---

*Built with passion for quantum computing and bioinformatics.*
