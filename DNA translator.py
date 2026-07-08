#!/usr/bin/env python3
"""
Quantum DNA Translator – 11‑Qubit Simulator (Silent, Fully Working)
---------------------------------------------------------------------
• 16 genetic codes, translation, CRISPR, 3D protein, State Viewer, etc.
• Guide RNA typed directly in the GUI – no pop‑up issues
• Zero‑length Bloch vectors drawn as dots to avoid runtime warnings
• Professional dark‑theme layout with no overlapping widgets

Author: Dr. Muhammad Saddam Khokhar – Quantummind.AI
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.widgets import TextBox, Button, CheckButtons, RadioButtons
import time, os

# =============================================================================
# Genetic data (16 tables)
# =============================================================================
_AMINO = ["A","C","D","E","F","G","H","I","K","L",
          "M","N","P","Q","R","S","T","V","W","Y","*"]
_AA2IDX = {aa:i for i,aa in enumerate(_AMINO)}
_AA2IDX_INV = {v:k for k,v in enumerate(_AMINO)}

_THREE_LETTER = {
    'A':'Ala','C':'Cys','D':'Asp','E':'Glu','F':'Phe','G':'Gly','H':'His',
    'I':'Ile','K':'Lys','L':'Leu','M':'Met','N':'Asn','P':'Pro','Q':'Gln',
    'R':'Arg','S':'Ser','T':'Thr','V':'Val','W':'Trp','Y':'Tyr','*':'Stop'
}

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

_COMPLEMENT = {'A':'T','T':'A','G':'C','C':'G'}
def reverse_complement(seq):
    return ''.join(_COMPLEMENT.get(base, base) for base in reversed(seq))

def codon_to_bits(codon):
    m = {'A':0b00, 'T':0b01, 'G':0b10, 'C':0b11}
    bits = 0
    for n in codon:
        bits = (bits << 2) | m[n]
    return bits

def build_unitary_from_code(code_name):
    if code_name not in GENETIC_CODES:
        code_name = "Standard"
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

class QuantumState:
    def __init__(self, n=11, noise_prob=0.0):
        self.n = n
        self.dim = 1 << n
        self.state = np.zeros(self.dim, dtype=complex)
        self.state[0] = 1.0
        self.noise_prob = noise_prob

    def apply_depolarizing_noise(self, qubit):
        if self.noise_prob <= 0 or np.random.rand() > self.noise_prob:
            return
        r = np.random.rand()
        if r < self.noise_prob/3:
            X = np.array([[0,1],[1,0]], dtype=complex)
            self.apply_gate(X, [qubit])
        elif r < 2*self.noise_prob/3:
            Y = np.array([[0,-1j],[1j,0]], dtype=complex)
            self.apply_gate(Y, [qubit])
        elif r < self.noise_prob:
            Z = np.array([[1,0],[0,-1]], dtype=complex)
            self.apply_gate(Z, [qubit])

    def apply_gate(self, gate, qubits):
        k = len(qubits)
        gate_t = gate.reshape([2]*(2*k))
        tensor = self.state.reshape([2]*self.n)
        others = [q for q in range(self.n) if q not in qubits]
        perm = list(qubits) + others
        tensor = np.transpose(tensor, perm)
        tensor = np.tensordot(gate_t, tensor, axes=[list(range(k,2*k)), list(range(k))])
        inv = np.argsort(perm)
        tensor = np.transpose(tensor, inv)
        self.state = tensor.reshape(-1)
        for q in qubits:
            self.apply_depolarizing_noise(q)

    def set_qubit(self, q, val):
        if val:
            X = np.array([[0,1],[1,0]], dtype=complex)
            self.apply_gate(X, [q])

    def get_probabilities(self, qubits):
        k = len(qubits)
        probs = np.zeros(1<<k)
        full = np.abs(self.state)**2
        for idx, amp in enumerate(full):
            bits = 0
            for i, q in enumerate(qubits):
                if (idx >> q) & 1:
                    bits |= (1<<i)
            probs[bits] += amp
        return probs

    def measure(self, qubits):
        probs = self.get_probabilities(qubits)
        return np.random.choice(len(probs), p=probs)

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
        sx = np.array([[0,1],[1,0]], dtype=complex)
        sy = np.array([[0,-1j],[1j,0]], dtype=complex)
        sz = np.array([[1,0],[0,-1]], dtype=complex)
        x = np.real(np.trace(rho @ sx))
        y = np.real(np.trace(rho @ sy))
        z = np.real(np.trace(rho @ sz))
        return x, y, z

# =============================================================================
# Main GUI
# =============================================================================
class QuantumDNAGui:
    def __init__(self):
        self.bg_main = '#0B0C10'
        self.bg_panel = '#13161D'
        self.cyan = '#00FFFF'
        self.orange = '#FF4500'
        self.gold = '#FFD700'
        self.white = '#FFFFFF'

        self.fig = plt.figure(figsize=(28, 20), facecolor=self.bg_main)
        self.fig.patch.set_facecolor(self.bg_main)

        outer = GridSpec(2, 2, figure=self.fig, left=0.04, right=0.98,
                         top=0.93, bottom=0.05, width_ratios=[1.8, 2.2],
                         height_ratios=[1,1], wspace=0.25, hspace=0.3)

        self.fig.text(0.5, 0.97, 'QUANTUM DNA TRANSLATOR – Dr. Muhammad Saddam Khokhar',
                      ha='center', va='top', fontsize=18, color=self.cyan, fontweight='bold')

        # Left column: circuit and probabilities
        self.ax_circuit = self.fig.add_subplot(outer[0,0], facecolor=self.bg_panel)
        self.ax_circuit.set_title('TRANSLATION CIRCUIT', color=self.cyan, fontsize=14, pad=10)
        self.ax_circuit.set_xlim(0,10); self.ax_circuit.set_ylim(-1,12)
        self.ax_circuit.axis('off')

        self.ax_prob = self.fig.add_subplot(outer[1,0], facecolor=self.bg_panel)
        self.ax_prob.set_title('OUTPUT QUBITS MEASUREMENT PROBABILITIES',
                               color=self.cyan, fontsize=12, pad=8)
        self.prob_bars = self.ax_prob.bar(range(32), [0]*32, color=self.orange,
                                          edgecolor='white', alpha=0.8)
        self.ax_prob.set_xlim(-0.5,31.5); self.ax_prob.set_ylim(0,1)
        self.ax_prob.set_xlabel('Amino Acid Index (0-20, others unused)', color=self.white, fontsize=10)
        self.ax_prob.set_ylabel('Probability', color=self.white, fontsize=10)
        self.ax_prob.tick_params(colors=self.white, labelsize=8)
        self.ax_prob.set_facecolor(self.bg_panel)

        # Right column: Bloch + controls
        right = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[:,1],
                                        height_ratios=[1.2, 2.2], hspace=0.3)

        # Bloch spheres
        bloch_gs = GridSpecFromSubplotSpec(2, 3, subplot_spec=right[0],
                                           wspace=0.4, hspace=0.5)
        self.bloch_axes = []
        for r in range(2):
            for c in range(3):
                ax = self.fig.add_subplot(bloch_gs[r,c], projection='3d',
                                          facecolor=self.bg_panel)
                ax.set_title(f'Qubit {r*3+c}', color=self.cyan, fontsize=10)
                ax.set_xlim(-1.2,1.2); ax.set_ylim(-1.2,1.2); ax.set_zlim(-1.2,1.2)
                ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
                ax.grid(False)
                u = np.linspace(0, 2*np.pi, 20)
                v = np.linspace(0, np.pi, 20)
                x = np.outer(np.cos(u), np.sin(v))
                y = np.outer(np.sin(u), np.sin(v))
                z = np.outer(np.ones(np.size(u)), np.cos(v))
                ax.plot_wireframe(x,y,z, color='grey', alpha=0.3, linewidth=0.3)
                ax.quiver(0,0,0, 1.5,0,0, color='red', arrow_length_ratio=0.1)
                ax.quiver(0,0,0, 0,1.5,0, color='green', arrow_length_ratio=0.1)
                ax.quiver(0,0,0, 0,0,1.5, color='blue', arrow_length_ratio=0.1)
                # initial Bloch vector (|0⟩)
                ax.quiver(0,0,0, 0,0,1, color=self.gold, lw=2, arrow_length_ratio=0.15)
                self.bloch_axes.append(ax)

        # Controls – 8 rows (perfectly spaced)
        ctrl_gs = GridSpecFromSubplotSpec(8, 1, subplot_spec=right[1],
                                          height_ratios=[0.03, 0.06, 0.04, 0.06, 0.22, 0.06, 0.22, 0.08],
                                          hspace=0.12)

        # Row0: DNA input label
        ax_lbl = self.fig.add_subplot(ctrl_gs[0], facecolor=self.bg_panel)
        ax_lbl.axis('off')
        ax_lbl.text(0.5, 0.5, 'Enter DNA sequence:', color=self.white, fontsize=11, ha='center', va='center')

        # Row1: DNA TextBox
        self.ax_text = self.fig.add_subplot(ctrl_gs[1], facecolor=self.bg_panel)
        self.textbox = TextBox(self.ax_text, '', initial='ATGTTTGGCTAATGCCGTACGTAGCTAGCTAGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG')
        self.textbox.on_submit(self.run_translation)

        # Row2: CRISPR guide RNA + Find button
        ax_crispr = self.fig.add_subplot(ctrl_gs[2], facecolor=self.bg_panel)
        ax_crispr.axis('off')
        guide_ax = ax_crispr.inset_axes([0.02, 0.15, 0.55, 0.7])
        self.guide_textbox = TextBox(guide_ax, 'Guide RNA (20 nt):', initial='TAGCTAGCTAGCATGCATGC')
        self.guide_textbox.label.set_color(self.white)
        self.guide_textbox.label.set_fontsize(8)
        find_ax = ax_crispr.inset_axes([0.60, 0.15, 0.2, 0.7])
        self.btn_find = Button(find_ax, 'Find Off-Targets',
                               color='#2A2F3A', hovercolor='#3E4450')
        self.btn_find.label.set_color(self.cyan); self.btn_find.label.set_fontsize(7)
        self.btn_find.on_clicked(self.run_crispr_search)

        # Row3: Translate button + timing
        ax_btn = self.fig.add_subplot(ctrl_gs[3], facecolor=self.bg_panel)
        ax_btn.axis('off')
        self.btn_ax = ax_btn.inset_axes([0.1, 0.1, 0.4, 0.8])
        self.btn_translate = Button(self.btn_ax, 'Translate',
                                    color='#2A2F3A', hovercolor='#3E4450')
        self.btn_translate.label.set_color(self.cyan)
        self.btn_translate.label.set_fontsize(12)
        self.btn_translate.on_clicked(self.run_translation)
        self.timing_label = ax_btn.text(0.55, 0.5, 'Time: --', color=self.gold, fontsize=10, va='center')

        # Row4: Protein / Result
        ax_prot = self.fig.add_subplot(ctrl_gs[4], facecolor=self.bg_panel)
        ax_prot.set_title('PROTEIN / RESULT', color=self.cyan, fontsize=12, pad=2)
        self.protein_text = ax_prot.text(0.05, 0.7, '', color=self.gold,
                                         fontsize=14, fontfamily='monospace', va='center')
        self.remark_text = ax_prot.text(0.05, 0.4, '', color=self.orange, fontsize=10, va='center')
        ax_prot.axis('off')

        # Row5: Combined row for Display Opts, Report, Noise
        ax_combined = self.fig.add_subplot(ctrl_gs[5], facecolor=self.bg_panel)
        ax_combined.axis('off')
        disp_ax = ax_combined.inset_axes([0.02, 0.1, 0.22, 0.8])
        self.btn_display = Button(disp_ax, 'Display Opts',
                                  color='#2A2F3A', hovercolor='#3E4450')
        self.btn_display.label.set_color(self.cyan); self.btn_display.label.set_fontsize(9)
        self.btn_display.on_clicked(self.show_display_options)
        rep_ax = ax_combined.inset_axes([0.28, 0.1, 0.22, 0.8])
        self.btn_report = Button(rep_ax, 'Report', color='#1E2230', hovercolor='#3E4450')
        self.btn_report.label.set_color(self.gold); self.btn_report.label.set_fontsize(11)
        self.btn_report.on_clicked(self.generate_report)
        noise_ax = ax_combined.inset_axes([0.55, 0.1, 0.40, 0.8])
        self.chk_noise = CheckButtons(noise_ax, ['Depol. Noise (p=0.01)'], [False])
        self.chk_noise.labels[0].set_color(self.white); self.chk_noise.labels[0].set_fontsize(10)
        self.noise_enabled = False
        def toggle_noise(label):
            self.noise_enabled = not self.noise_enabled
        self.chk_noise.on_clicked(toggle_noise)

        # Row6: Genetic code selection – 4×4 grid
        ax_gc = self.fig.add_subplot(ctrl_gs[6], facecolor=self.bg_panel)
        ax_gc.set_title('Select Genetic Code', color=self.cyan, fontsize=12, pad=4)
        ax_gc.axis('off')
        code_names = list(GENETIC_CODES.keys())
        short_names = [
            "Standard", "Vertebrate mito", "Yeast mito", "Mold/Protoz mito",
            "Invert mito", "Ciliate/Hexamita", "Echinoderm mito", "Euplotid nuclear",
            "Alt yeast", "Ascidian mito", "Alt flatworm", "Blepharisma",
            "Chlorophycean", "Trematode mito", "Scenedesmus", "Pterobranchia"
        ]
        n_cols = 4
        btn_w = 0.21; btn_h = 0.16
        start_x = 0.02; start_y = 0.78
        gap_x = 0.03; gap_y = 0.10
        self.code_buttons = {}
        for idx, code_name in enumerate(code_names):
            r = idx // n_cols; c = idx % n_cols
            left = start_x + c*(btn_w + gap_x)
            bottom = start_y - r*(btn_h + gap_y)
            ax_btn = ax_gc.inset_axes([left, bottom, btn_w, btn_h])
            btn = Button(ax_btn, short_names[idx], color='#2A2F3A', hovercolor='#3E4450')
            btn.label.set_color(self.cyan); btn.label.set_fontsize(9)
            btn.on_clicked(lambda event, code=code_name: self.set_genetic_code(code))
            self.code_buttons[code_name] = btn
        self.gc_label = ax_gc.text(0.5, 0.04, 'Current: Standard', color=self.gold,
                                   fontsize=10, ha='center', va='center', transform=ax_gc.transAxes)

        # Row7: Small buttons (State, Heatmap, DNA Img, 3D Plot)
        ax_row7 = self.fig.add_subplot(ctrl_gs[7], facecolor=self.bg_panel)
        ax_row7.axis('off')
        positions = [0.05, 0.28, 0.51, 0.74]
        labels = ['State', 'Heatmap', 'DNA Img', '3D Plot']
        callbacks = [self.show_state_viewer, self.show_genetic_code_viz,
                     self.show_dna_image, self.show_protein_3d]
        self.control_buttons = {}
        for i in range(4):
            ax_b = ax_row7.inset_axes([positions[i], 0.1, 0.18, 0.8])
            btn = Button(ax_b, labels[i], color='#2A2F3A', hovercolor='#3E4450')
            btn.label.set_color(self.cyan); btn.label.set_fontsize(10)
            btn.on_clicked(callbacks[i])
            self.control_buttons[labels[i]] = btn

        # State variables
        self.codons = []
        self.protein = ""
        self.translation_done = False
        self.start_time = None
        self.noise_prob = 0.01
        self.current_code = "Standard"
        self.U_TRANS = build_unitary_from_code(self.current_code)
        self.display_opts = {'verbose': False, 'spaces': True, 'show_dna': False, 'strand': 'forward'}

    # ---------- Helper methods ----------
    def set_strand(self, label):
        self.display_opts['strand'] = label

    def set_genetic_code(self, name):
        if name in GENETIC_CODES:
            self.current_code = name
            self.U_TRANS = build_unitary_from_code(name)
            self.gc_label.set_text(f'Current: {name}')
            self.fig.canvas.draw_idle()

    def draw_circuit(self, bits):
        self.ax_circuit.clear()
        self.ax_circuit.set_xlim(0,10); self.ax_circuit.set_ylim(-1,12)
        self.ax_circuit.axis('off')
        for i in range(6):
            y = 10 - i
            self.ax_circuit.plot([1,3.5], [y,y], color='white', lw=1)
            self.ax_circuit.text(0.5, y, f'q{i}', color=self.white, fontsize=8)
            bit = (bits >> (5-i)) & 1
            self.ax_circuit.text(2.2, y+0.3, f'|{bit}⟩', color=self.cyan, fontsize=7)
        for i in range(5):
            y = 4 - i
            self.ax_circuit.plot([1,3.5], [y,y], color='white', lw=1)
            self.ax_circuit.text(0.5, y, f'q{6+i}', color=self.white, fontsize=8)
            self.ax_circuit.text(2.2, y+0.3, '|0⟩', color=self.cyan, fontsize=7)
        rect = plt.Rectangle((3.5,3), 1.5, 8, fill=False, edgecolor=self.gold, lw=2)
        self.ax_circuit.add_patch(rect)
        self.ax_circuit.text(4.25, 7, 'U', color=self.gold, fontsize=14, ha='center', va='center')
        for i in range(5):
            y = 4 - i
            self.ax_circuit.plot([5.5,6.5], [y,y], color='white', lw=1)
            self.ax_circuit.add_patch(plt.Circle((7.0, y), 0.3, fill=False, edgecolor=self.orange))
            self.ax_circuit.plot([7.5,8.5], [y,y], color='white', lw=1)
            self.ax_circuit.text(9.0, y, f'a{i}', color=self.white, fontsize=8)
        self.fig.canvas.draw_idle()

    def update_bloch(self, state):
        """Update Bloch spheres – draws a dot for zero‑length vectors to avoid warnings."""
        for q in range(6):
            ax = self.bloch_axes[q]
            x, y, z = state.bloch_vector(q)
            for coll in ax.collections:
                coll.remove()
            # Remove previous quiver (but keep fixed axis arrows – they are not collections)
            # To remove old golden arrow, we must clear all children? We only remove collections
            # which are quivers. The scatter points are also collections, so they get removed too.
            if np.linalg.norm([x, y, z]) < 1e-9:
                ax.scatter(0, 0, 0, color=self.gold, s=30)
            else:
                ax.quiver(0, 0, 0, x, y, z, color=self.gold, lw=2, arrow_length_ratio=0.15)
        self.fig.canvas.draw_idle()

    def update_probs(self, probs):
        for bar, h in zip(self.prob_bars, probs):
            bar.set_height(h)
        self.fig.canvas.draw_idle()

    # ---------- Translation ----------
    def process_codon(self, codon):
        qs = QuantumState(11, noise_prob=self.noise_prob if self.noise_enabled else 0.0)
        bits = codon_to_bits(codon)
        for i in range(6):
            qs.set_qubit(i, (bits >> (5-i)) & 1)
        self.draw_circuit(bits)
        self.update_bloch(qs)
        plt.pause(0.3)
        qs.apply_gate(self.U_TRANS, list(range(11)))
        self.update_bloch(qs)
        plt.pause(0.3)
        probs = qs.get_probabilities([6,7,8,9,10])
        self.update_probs(probs)
        outcome = qs.measure([6,7,8,9,10])
        idx = outcome if outcome < 21 else 20
        amino = _AMINO[idx]
        return amino

    def _format_amino(self, aa):
        if self.display_opts['verbose']:
            return _THREE_LETTER.get(aa, aa)
        return aa

    def run_translation(self, event=None):
        if self.translation_done:
            self.reset_translation()
            return
        dna = self.textbox.text.strip().upper()
        dna = ''.join(c for c in dna if c in 'ATGC')
        if len(dna) % 3 != 0:
            dna = dna[:-(len(dna)%3)]
        if not dna:
            return
        strand = self.display_opts['strand']
        if strand == 'reverse':
            seq = reverse_complement(dna)
        elif strand == 'both':
            seq = dna
        else:
            seq = dna
        self.codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
        self.protein = ""
        self.protein_text.set_text("")
        self.remark_text.set_text("")
        self.start_time = time.time()
        for codon in self.codons:
            aa = self.process_codon(codon)
            if aa == '*':
                self.protein += '*'
                self.remark_text.set_text('Stop codon encountered.')
                break
            disp_aa = self._format_amino(aa)
            if self.display_opts['spaces']:
                disp_aa += ' '
            self.protein += disp_aa
            if self.display_opts['show_dna']:
                display = f"DNA: {seq}\nProtein: {self.protein}"
            else:
                display = self.protein
            self.protein_text.set_text(display)
            self.fig.canvas.draw_idle()
            time.sleep(0.02)
        else:
            self.remark_text.set_text('Translation complete.')
        elapsed = time.time() - self.start_time
        self.timing_label.set_text(f'Time: {elapsed:.2f}s')
        self.btn_translate.label.set_text('New Translation')
        self.btn_translate.on_clicked(self.run_translation)
        self.translation_done = True

    def reset_translation(self):
        self.protein = ""
        self.protein_text.set_text("")
        self.remark_text.set_text("")
        self.codons = []
        self.translation_done = False
        self.start_time = None
        self.timing_label.set_text('Time: --')
        self.ax_circuit.clear()
        self.ax_circuit.set_xlim(0,10); self.ax_circuit.set_ylim(-1,12)
        self.ax_circuit.axis('off')
        for bar in self.prob_bars:
            bar.set_height(0)
        for ax in self.bloch_axes:
            for coll in ax.collections:
                coll.remove()
            ax.quiver(0,0,0, 0,0,1, color=self.gold, lw=2, arrow_length_ratio=0.15)
        self.btn_translate.label.set_text('Translate')
        self.btn_translate.on_clicked(self.run_translation)
        self.fig.canvas.draw_idle()

    # ---------- CRISPR search ----------
    def run_crispr_search(self, event=None):
        guide = self.guide_textbox.text.strip().upper()[:20]
        if len(guide) < 20:
            self.remark_text.set_text('Guide RNA must be 20 nt.')
            return
        target_dna = self.textbox.text.strip().upper()
        target_dna = ''.join(c for c in target_dna if c in 'ATGC')
        mismatches = []
        for i in range(len(target_dna) - 19):
            window = target_dna[i:i+20]
            dist = sum(1 for a, b in zip(guide, window) if a != b)
            if dist <= 3:
                mismatches.append((i, dist, window))
        fig6, ax = plt.subplots(figsize=(10, 6), facecolor='#0B0C10')
        fig6.patch.set_facecolor('#0B0C10')
        ax.set_facecolor('#0B0C10')
        ax.axis('off')
        text = f'CRISPR Off‑Target Results\n\nTarget: {len(target_dna)} bp\nGuide: {guide}\n\n'
        text += f'Off‑target sites (≤3 mismatches): {len(mismatches)} found\n'
        if mismatches:
            for pos, dist, seq in mismatches[:20]:
                text += f'  Pos {pos}: {seq} ({dist} mm)\n'
        else:
            text += '  None found.\n'
        ax.text(0.05, 0.95, text, color=self.white, fontsize=10, va='top', fontfamily='monospace')
        plt.show(block=False)

    # ---------- Report generation ----------
    def generate_report(self, event=None):
        dna = self.textbox.text.strip().upper()
        protein = self.protein_text.get_text()
        remarks = self.remark_text.get_text()
        noise = self.noise_enabled
        elapsed = time.time() - self.start_time if self.start_time else 0
        report = f"""
Quantum DNA Translation Report
Generated by Dr. Muhammad Saddam Khokhar
-----------------------------------------
Date: {time.ctime()}
Input DNA: {dna}
Protein Sequence: {protein}
Remarks: {remarks}
Noise Model: {'Enabled (p=' + str(self.noise_prob) + ')' if noise else 'Disabled'}
Translation Time: {elapsed:.2f} seconds
Number of codons: {len(self.codons)}
Genetic Code: {self.current_code}
Display: {'Verbose' if self.display_opts['verbose'] else 'Compact'}, Spaces={'Yes' if self.display_opts['spaces'] else 'No'}, Strand={self.display_opts['strand']}
-----------------------------------------
"""
        with open('quantum_translation_report.txt', 'w') as f:
            f.write(report)
        self.show_popup_report(report)
        self.btn_report.label.set_text('Report Saved')
        self.fig.canvas.draw_idle()
        time.sleep(1)
        self.btn_report.label.set_text('Report')

    def show_popup_report(self, text):
        fig2, ax = plt.subplots(figsize=(8,6))
        fig2.patch.set_facecolor('#0B0C10')
        ax.set_facecolor('#0B0C10')
        ax.axis('off')
        ax.text(0.05, 0.9, text, fontfamily='monospace', color=self.white,
                va='top', fontsize=10)
        close_ax = plt.axes([0.8, 0.05, 0.15, 0.06])
        close_btn = Button(close_ax, 'Close', color='#2A2F3A', hovercolor='#3E4450')
        close_btn.label.set_color(self.cyan)
        def close(event):
            plt.close(fig2)
        close_btn.on_clicked(close)
        plt.show(block=False)

    # ---------- Display Options Popup ----------
    def show_display_options(self, event=None):
        fig4, axs = plt.subplots(3, 1, figsize=(8, 6))
        fig4.patch.set_facecolor('#0B0C10')
        for ax in axs:
            ax.set_facecolor('#0B0C10')
            ax.tick_params(colors='white')
        ax_verbose = axs[0]; ax_verbose.set_title('Residue format', color=self.white)
        chk_verbose = CheckButtons(ax_verbose, ['Verbose (three-letter)'], [self.display_opts['verbose']])
        chk_verbose.labels[0].set_color(self.white); chk_verbose.labels[0].set_fontsize(12)
        ax_spaces = axs[1]; ax_spaces.set_title('Spacing', color=self.white)
        chk_spaces = CheckButtons(ax_spaces, ['Spaces between residues'], [self.display_opts['spaces']])
        chk_spaces.labels[0].set_color(self.white); chk_spaces.labels[0].set_fontsize(12)
        ax_showdna = axs[2]; ax_showdna.set_title('Include DNA', color=self.white)
        chk_showdna = CheckButtons(ax_showdna, ['Show nucleotide sequence'], [self.display_opts['show_dna']])
        chk_showdna.labels[0].set_color(self.white); chk_showdna.labels[0].set_fontsize(12)
        strand_ax = plt.axes([0.1, 0.05, 0.8, 0.1])
        strand_ax.set_facecolor('#0B0C10')
        strand_radio = RadioButtons(strand_ax, ['forward', 'reverse', 'both'], active=0)
        for label in strand_radio.labels:
            label.set_color(self.white); label.set_fontsize(12)
        chosen_strand = [self.display_opts['strand']]
        def set_strand(label):
            chosen_strand[0] = label
        strand_radio.on_clicked(set_strand)
        apply_ax = plt.axes([0.7, 0.02, 0.2, 0.05])
        apply_btn = Button(apply_ax, 'Apply', color='#2A2F3A', hovercolor='#3E4450')
        apply_btn.label.set_color(self.cyan); apply_btn.label.set_fontsize(12)
        def apply_opts(event):
            self.display_opts['verbose'] = chk_verbose.get_status()[0]
            self.display_opts['spaces'] = chk_spaces.get_status()[0]
            self.display_opts['show_dna'] = chk_showdna.get_status()[0]
            self.display_opts['strand'] = chosen_strand[0]
            plt.close(fig4)
        apply_btn.on_clicked(apply_opts)
        plt.show(block=False)

    # ---------- State Viewer ----------
    def show_state_viewer(self, event=None):
        if not self.codons:
            return
        codon = self.codons[0]
        qs = QuantumState(11, noise_prob=self.noise_prob if self.noise_enabled else 0.0)
        bits = codon_to_bits(codon)
        for i in range(6):
            qs.set_qubit(i, (bits >> (5-i)) & 1)
        qs.apply_gate(self.U_TRANS, list(range(11)))
        fig2, axs = plt.subplots(2,2, figsize=(12,8))
        fig2.suptitle('Quantum State Viewer – First Codon', color='cyan')
        codon_indices = [codon_to_bits(c) for c in GENETIC_CODES[self.current_code]]
        amps = np.array([np.abs(qs.state[(i << 5)]) for i in codon_indices])
        axs[0,0].bar(range(64), amps, color='orange')
        axs[0,0].set_title('Codon State Amplitudes (64)', color='white')
        axs[0,0].set_xlabel('Codon index'); axs[0,0].set_ylabel('|Amplitude|')
        axs[0,0].set_xticks(range(0,64,8))
        rho0 = qs.get_reduced_density_matrix(0)
        axs[0,1].text(0.5, 0.5, f'Reduced density matrix\nQubit 0:\n\n{np.array2string(rho0, precision=2)}',
                      ha='center', va='center', color='white', fontsize=10,
                      fontfamily='monospace')
        axs[0,1].set_title('Qubit 0 Density', color='white')
        axs[1,0].set_title('Bloch Vectors', color='white')
        for q in range(6):
            x,y,z = qs.bloch_vector(q)
            if np.linalg.norm([x,y,z]) < 1e-9:
                axs[1,0].scatter(0,0, color='cyan')
            else:
                axs[1,0].quiver(0,0,x,y, color='cyan', label=f'Q{q}')
        axs[1,0].legend(); axs[1,0].set_xlim(-1.2,1.2); axs[1,0].set_ylim(-1.2,1.2)
        probs = qs.get_probabilities([6,7,8,9,10])
        axs[1,1].bar(range(32), probs, color='orange')
        axs[1,1].set_title('Output Probs', color='white')
        for ax in axs.flat:
            ax.set_facecolor('#0B0C10')
            ax.tick_params(colors='white')
        fig2.patch.set_facecolor('#0B0C10')
        plt.show(block=False)

    # ---------- Genetic code heatmap ----------
    def show_genetic_code_viz(self, event=None):
        bases = ['T','C','A','G']
        fig3, ax = plt.subplots(figsize=(8,6))
        fig3.patch.set_facecolor('#0B0C10')
        ax.set_facecolor('#0B0C10')
        ax.set_title('Genetic Code Heatmap', color='cyan', fontsize=14)
        grid = [['' for _ in range(4)] for _ in range(4)]
        for i, row_base in enumerate(bases):
            for j, col_base in enumerate(bases):
                codon = row_base + col_base + 'T'
                aa = GENETIC_CODES[self.current_code].get(codon, '?')
                grid[i][j] = aa
        colors = plt.cm.tab20(np.linspace(0,1,20))
        for i in range(4):
            for j in range(4):
                aa = grid[i][j]
                idx = _AA2IDX.get(aa, 0)
                ax.add_patch(plt.Rectangle((j, 3-i), 1, 1, color=colors[idx%20]))
                ax.text(j+0.5, 3-i+0.5, aa, ha='center', va='center', fontsize=12)
        ax.set_xlim(0,4); ax.set_ylim(0,4)
        ax.set_xticks([0.5,1.5,2.5,3.5]); ax.set_xticklabels(bases, color='white')
        ax.set_yticks([0.5,1.5,2.5,3.5]); ax.set_yticklabels(bases[::-1], color='white')
        ax.tick_params(colors='white')
        plt.show(block=False)

    # ---------- DNA Image Generator ----------
    def show_dna_image(self, event=None):
        dna = self.textbox.text.strip().upper()
        dna = ''.join(c for c in dna if c in 'ATGC')
        if not dna:
            return
        fig5, ax = plt.subplots(figsize=(12, 2))
        fig5.patch.set_facecolor('#0B0C10')
        ax.set_facecolor('#0B0C10')
        colors = {'A': 'red', 'T': 'blue', 'G': 'green', 'C': 'yellow'}
        for i, base in enumerate(dna):
            ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=colors[base]))
            ax.text(i+0.5, 0.5, base, ha='center', va='center', fontweight='bold')
        ax.set_xlim(0, len(dna))
        ax.set_ylim(0, 1)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title('DNA Sequence Visualization', color=self.cyan, fontsize=14)
        plt.show(block=False)

    # ---------- 3D Protein Structure Placeholder ----------
    def show_protein_3d(self, event=None):
        protein = self.protein_text.get_text()
        residues = [aa for aa in protein if aa != ' ']
        if not residues:
            return
        fig7 = plt.figure(figsize=(10,8))
        ax = fig7.add_subplot(111, projection='3d', facecolor='#0B0C10')
        fig7.patch.set_facecolor('#0B0C10')
        ax.set_title('Protein 3D Structure (Random Coil)', color=self.cyan, fontsize=14)
        n = len(residues)
        phi = np.linspace(0, 4*np.pi, n)
        r = np.linspace(0.5, 2.0, n)
        x = r * np.cos(phi) + np.random.normal(0, 0.1, n)
        y = r * np.sin(phi) + np.random.normal(0, 0.1, n)
        z = np.linspace(0, 2.5, n) + np.random.normal(0, 0.1, n)
        color_map = {'A':'red','C':'yellow','D':'blue','E':'blue','F':'green',
                     'G':'red','H':'blue','I':'green','K':'blue','L':'green',
                     'M':'green','N':'blue','P':'red','Q':'blue','R':'blue',
                     'S':'red','T':'red','V':'green','W':'green','Y':'green',
                     '*':'white'}
        colors = [color_map.get(aa, 'grey') for aa in residues]
        ax.scatter(x, y, z, c=colors, s=50)
        for i, (xi,yi,zi,aa) in enumerate(zip(x,y,z,residues)):
            ax.text(xi, yi, zi, aa, color='white', fontsize=8)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
        ax.xaxis.label.set_color('white'); ax.yaxis.label.set_color('white'); ax.zaxis.label.set_color('white')
        ax.tick_params(colors='white')
        plt.show(block=False)

if __name__ == '__main__':
    gui = QuantumDNAGui()
    plt.show()