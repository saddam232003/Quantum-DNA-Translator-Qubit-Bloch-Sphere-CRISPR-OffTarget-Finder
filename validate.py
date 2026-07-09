#!/usr/bin/env python3
"""
Headless validation & benchmarking of the Quantum DNA Translator core logic.
Extracts the exact numerical routines from the submitted source (genetic code
tables, codon_to_bits, build_unitary_from_code, QuantumState) with all
matplotlib/GUI code stripped out, and subjects them to four independent
checks:
  1. Unitarity of the constructed "translation operator" U, for all 16
     genetic code tables, for both the ORIGINAL construction routine and a
     CORRECTED construction routine.
  2. Round-trip correctness of codon -> amino-acid translation under U,
     cross-checked against Biopython's independent codon table implementation.
  3. Bloch-vector / reduced-density-matrix sanity checks on known states.
  4. Wall-clock benchmarking of the quantum-simulated translation pathway
     against a plain classical dictionary lookup, as a function of sequence
     length.
"""
import json
import random
import time

import numpy as np
from Bio.Data import CodonTable

# ----------------------------------------------------------------------
# 1. Exact genetic-code data copied verbatim from the submitted source
# ----------------------------------------------------------------------
_AMINO = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
          "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "*"]
_AA2IDX = {aa: i for i, aa in enumerate(_AMINO)}

GENETIC_CODES = {
    "Standard": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Vertebrate mitochondrial": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'M','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'*','AGG':'*',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Yeast mitochondrial": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'T','CTC':'T','CTA':'T','CTG':'T','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'M','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Mold/Protozoan/Coelenterate mito": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Invertebrate mitochondrial": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'M','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'S','AGG':'S',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Ciliate/Dasycladacean/Hexamita nuclear": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'Q','TAG':'Q','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Echinoderm/Flatworm mitochondrial": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'N','AAG':'K','AGT':'S','AGC':'S','AGA':'S','AGG':'S',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Euplotid nuclear": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'C','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Alternative yeast nuclear": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Ascidian mitochondrial": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'M','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'G','AGG':'G',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Alternative flatworm mitochondrial": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'N','AAG':'K','AGT':'S','AGC':'S','AGA':'S','AGG':'S',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Blepharisma nuclear": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Chlorophycean mitochondrial": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Trematode mitochondrial": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'M','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'N','AAG':'K','AGT':'S','AGC':'S','AGA':'S','AGG':'S',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Scenedesmus obliquus mito": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    },
    "Pterobranchia mitochondrial": {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'S','AGG':'K',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
    }
}

# Map our table names onto Biopython's NCBI transl_table IDs for cross-checking
NCBI_ID_MAP = {
    "Standard": 1, "Vertebrate mitochondrial": 2, "Yeast mitochondrial": 3,
    "Mold/Protozoan/Coelenterate mito": 4, "Invertebrate mitochondrial": 5,
    "Ciliate/Dasycladacean/Hexamita nuclear": 6, "Echinoderm/Flatworm mitochondrial": 9,
    "Euplotid nuclear": 10, "Alternative yeast nuclear": 12, "Ascidian mitochondrial": 13,
    "Alternative flatworm mitochondrial": 14, "Blepharisma nuclear": 15,
    "Chlorophycean mitochondrial": 16, "Trematode mitochondrial": 21,
    "Scenedesmus obliquus mito": 22, "Pterobranchia mitochondrial": 24,
}


def codon_to_bits(codon):
    m = {'A': 0b00, 'T': 0b01, 'G': 0b10, 'C': 0b11}
    bits = 0
    for n in codon:
        bits = (bits << 2) | m[n]
    return bits


def build_unitary_from_code(code_name):
    """Verbatim port of the ORIGINAL construction routine from the source file."""
    codon_table = GENETIC_CODES[code_name]
    dim = 1 << 11
    U = np.zeros((dim, dim), dtype=complex)
    for codon, aa in codon_table.items():
        x = codon_to_bits(codon)
        y = _AA2IDX[aa] if aa != '*' else _AA2IDX['*']
        in_idx = (x << 5)
        out_idx = (x << 5) | y
        U[out_idx, in_idx] = 1.0
    for i in range(dim):
        if not np.any(U[:, i]):
            U[i, i] = 1.0
    return U


def build_unitary_from_code_fixed(code_name):
    """Corrected construction: allocate a bijection from unused columns to
    unused rows, instead of blindly setting U[i, i] = 1 for every column
    that has not yet been touched (which can collide with rows already
    claimed by a genuine codon mapping)."""
    codon_table = GENETIC_CODES[code_name]
    dim = 1 << 11
    U = np.zeros((dim, dim), dtype=complex)
    used_rows, assigned_cols = set(), set()
    for codon, aa in codon_table.items():
        x = codon_to_bits(codon)
        y = _AA2IDX[aa]
        in_idx, out_idx = (x << 5), (x << 5) | y
        U[out_idx, in_idx] = 1.0
        used_rows.add(out_idx)
        assigned_cols.add(in_idx)
    free_rows = sorted(set(range(dim)) - used_rows)
    remaining_cols = [i for i in range(dim) if i not in assigned_cols]
    assert len(free_rows) == len(remaining_cols)
    for col, row in zip(remaining_cols, free_rows):
        U[row, col] = 1.0
    return U


def check_unitary(U, label):
    dim = U.shape[0]
    dev = float(np.linalg.norm(U.conj().T @ U - np.eye(dim)))
    row_nnz = np.sum(np.abs(U) > 1e-9, axis=1)
    col_nnz = np.sum(np.abs(U) > 1e-9, axis=0)
    return {
        "label": label,
        "unitarity_frobenius_dev": dev,
        "rows_with_wrong_nnz": int(np.sum(row_nnz != 1)),
        "cols_with_wrong_nnz": int(np.sum(col_nnz != 1)),
        "max_row_weight": int(row_nnz.max()),
        "rows_with_zero_entries": int(np.sum(row_nnz == 0)),
    }


# ----------------------------------------------------------------------
# QuantumState -- verbatim core (matplotlib-dependent methods stripped)
# ----------------------------------------------------------------------
class QuantumState:
    def __init__(self, n=11):
        self.n = n
        self.dim = 1 << n
        self.state = np.zeros(self.dim, dtype=complex)
        self.state[0] = 1.0

    def apply_gate(self, gate, qubits):
        k = len(qubits)
        gate_t = gate.reshape([2] * (2 * k))
        tensor = self.state.reshape([2] * self.n)
        others = [q for q in range(self.n) if q not in qubits]
        perm = list(qubits) + others
        tensor = np.transpose(tensor, perm)
        tensor = np.tensordot(gate_t, tensor, axes=[list(range(k, 2 * k)), list(range(k))])
        inv = np.argsort(perm)
        tensor = np.transpose(tensor, inv)
        self.state = tensor.reshape(-1)

    def set_qubit(self, q, val):
        if val:
            X = np.array([[0, 1], [1, 0]], dtype=complex)
            self.apply_gate(X, [q])

    def get_reduced_density_matrix(self, qubit):
        n = self.n
        rho = np.outer(self.state, self.state.conj())
        rho = rho.reshape([2] * (2 * n))
        others = [i for i in range(n) if i != qubit]
        perm = [qubit] + others + [qubit + n] + [i + n for i in others]
        rho = np.transpose(rho, perm)
        env_dim = 1 << (n - 1)
        rho = rho.reshape((2, env_dim, 2, env_dim))
        rho = np.trace(rho, axis1=1, axis2=3)
        return rho

    def bloch_vector(self, qubit):
        rho = self.get_reduced_density_matrix(qubit)
        sx = np.array([[0, 1], [1, 0]], dtype=complex)
        sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
        sz = np.array([[1, 0], [0, -1]], dtype=complex)
        return (float(np.real(np.trace(rho @ sx))),
                float(np.real(np.trace(rho @ sy))),
                float(np.real(np.trace(rho @ sz))))


def main():
    report = {}

    # ---- Check 1 & 2: unitarity, original vs. fixed, across all 16 tables ----
    orig_results, fixed_results = [], []
    for name in GENETIC_CODES:
        orig_results.append(check_unitary(build_unitary_from_code(name), name))
        fixed_results.append(check_unitary(build_unitary_from_code_fixed(name), name))
    report["unitarity_original"] = orig_results
    report["unitarity_fixed"] = fixed_results
    report["n_tables_nonunitary_original"] = sum(
        1 for r in orig_results if r["unitarity_frobenius_dev"] > 1e-6)
    report["n_tables_unitary_fixed"] = sum(
        1 for r in fixed_results if r["unitarity_frobenius_dev"] < 1e-6)

    # ---- Check 3: translation correctness for the tool's ACTUAL usage
    #      pattern (basis-state input), cross-checked against Biopython ----
    total_checks, internal_failures, biopython_mismatches = 0, 0, 0
    for name, table in GENETIC_CODES.items():
        U = build_unitary_from_code(name)
        bio_table = CodonTable.unambiguous_dna_by_id[NCBI_ID_MAP[name]].forward_table
        bio_stops = set(CodonTable.unambiguous_dna_by_id[NCBI_ID_MAP[name]].stop_codons)
        for codon, aa in table.items():
            x = codon_to_bits(codon)
            vec_in = np.zeros(2048, dtype=complex)
            vec_in[x << 5] = 1.0
            vec_out = U @ vec_in
            out_idx = int(np.argmax(np.abs(vec_out)))
            predicted_y = out_idx & 0b11111
            expected_y = _AA2IDX[aa]
            total_checks += 1
            if predicted_y != expected_y or abs(vec_out[out_idx]) < 0.999:
                internal_failures += 1
            # cross-check against Biopython's independently maintained tables
            if aa == '*':
                if codon not in bio_stops:
                    biopython_mismatches += 1
            else:
                if bio_table.get(codon) != aa:
                    biopython_mismatches += 1
    report["translation_checks_total"] = total_checks
    report["translation_internal_failures"] = internal_failures
    report["translation_biopython_mismatches"] = biopython_mismatches

    # ---- Check 4: Bloch-vector / density-matrix sanity ----
    qs0 = QuantumState(11)
    bloch0 = qs0.bloch_vector(0)
    qs1 = QuantumState(11)
    qs1.set_qubit(0, 1)
    bloch1 = qs1.bloch_vector(0)
    H = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
    qsplus = QuantumState(11)
    qsplus.apply_gate(H, [0])
    blochplus = qsplus.bloch_vector(0)
    rho_test = qsplus.get_reduced_density_matrix(0)
    report["bloch_zero_state"] = bloch0
    report["bloch_one_state"] = bloch1
    report["bloch_plus_state"] = blochplus
    report["density_trace_real"] = float(np.real(np.trace(rho_test)))
    report["density_hermiticity_dev"] = float(np.linalg.norm(rho_test - rho_test.conj().T))

    # ---- Check 5: runtime benchmark, classical dict lookup vs. the
    #      per-codon quantum-simulated pathway (fresh QuantumState + full
    #      11-qubit apply_gate call per codon, exactly as in process_codon) ----
    def classical_translate(dna, table):
        out = []
        for i in range(0, len(dna) - 2, 3):
            aa = table[dna[i:i + 3]]
            if aa == '*':
                break
            out.append(aa)
        return ''.join(out)

    def quantum_translate(dna, U):
        out = []
        for i in range(0, len(dna) - 2, 3):
            codon = dna[i:i + 3]
            qs = QuantumState(11)
            bits = codon_to_bits(codon)
            for b in range(6):
                qs.set_qubit(b, (bits >> (5 - b)) & 1)
            qs.apply_gate(U, list(range(11)))
            idx = int(np.argmax(np.abs(qs.state) ** 2))
            aa = _AMINO[idx & 0b11111]
            if aa == '*':
                break
            out.append(aa)
        return ''.join(out)

    random.seed(42)
    std_table = GENETIC_CODES["Standard"]
    non_stop_codons = [c for c, aa in std_table.items() if aa != '*']
    U_std = build_unitary_from_code("Standard")

    timing = []
    for n_codons in [10, 50, 100, 300, 600]:
        dna = ''.join(random.choice(non_stop_codons) for _ in range(n_codons))
        reps = 3
        t0 = time.perf_counter()
        for _ in range(reps):
            p_classical = classical_translate(dna, std_table)
        t_classical = (time.perf_counter() - t0) / reps

        t0 = time.perf_counter()
        for _ in range(reps):
            p_quantum = quantum_translate(dna, U_std)
        t_quantum = (time.perf_counter() - t0) / reps

        assert p_classical == p_quantum, f"Mismatch at n={n_codons}"
        timing.append({
            "n_codons": n_codons,
            "classical_seconds": t_classical,
            "quantum_seconds": t_quantum,
            "slowdown_factor": t_quantum / t_classical if t_classical > 0 else None,
        })
    report["timing"] = timing

    print(json.dumps(report, indent=2))
    with open(r"C:\Users\MSKHOKHAR_PK\PycharmProjects\PythonProject1\validation_report.json", "w") as f:
        json.dump(report, f, indent=2)


if __name__ == "__main__":
    main()
