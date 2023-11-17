# SimpleQuantumCircuit

[![Build Status](https://github.com/malaydasat@gmail.com/SimpleQuantumCircuit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/malaydasat@gmail.com/SimpleQuantumCircuit.jl/actions/workflows/CI.yml?query=branch%3Amain)

A simple module to build quantum circuits that is also easy to learn. It uses Julia's native types. The module is designed to be easily extendable and also perform calculation at higher precisions.
This module provides a Unitary Evolution Matrix "U" for a sequence of gates. The U can be applied to any wavefunction which is represented as an array

## Usage

### Prepare a circuit
Use function circuit(number of qubits, Sequence of gates)
This plots the circuit and returns the unitary evolution matrix

```
julia> U = circuit(2, 
                   H([1]),
                   CX([1,2]) )
   ┌──────┐┌──────┐
1--┤1 H   ├┤1 CX  ├
   └──────┘└──────┘
           ┌──────┐
2----------┤2 CX  ├
           └──────┘
4×4 Matrix{Float64}:
 0.707107   0.707107  0.0        0.0
 0.0        0.0       0.707107  -0.707107
 0.0        0.0       0.707107   0.707107
 0.707107  -0.707107  0.0        0.0
```

### Going higher precision

To use an 512 bits arbitray precision number instead of 64 bit floating point

```
julia> setType(BigFloat)
BigFloat

julia> setprecision(512)
512

julia> U1 = circuit(2, 
                   H([1]),
                   CX([1,2]) )            
   ┌──────┐┌──────┐
1--┤1 H   ├┤1 CX  ├
   └──────┘└──────┘
           ┌──────┐
2----------┤2 CX  ├
           └──────┘
4×4 Matrix{BigFloat}:
 0.707107   0.707107  0.0        0.0
 0.0        0.0       0.707107  -0.707107
 0.0        0.0       0.707107   0.707107
 0.707107  -0.707107  0.0        0.0

```

### Defining a custom matrix for a gate in a circuit

matrix '[1 1; 1 -1]' is used as gate named "MY"

```
julia> U = circuit(2, 
                   H([1]),
                   gate([1 1; 1 -1], [2], "MY") )  
   ┌──────┐        
1--┤1 H   ├--------
   └──────┘        
           ┌──────┐
2----------┤1 MY  ├
           └──────┘
4×4 Matrix{BigFloat}:
 0.707107   0.707107   0.707107   0.707107
 0.707107  -0.707107   0.707107  -0.707107
 0.707107   0.707107  -0.707107  -0.707107
 0.707107  -0.707107  -0.707107   0.707107
```
### working with gates

There are many predefined gates as Ugate, Igate, X, Y, Z, H, RX, RY, RZ, swap, cnot, XX, YY, ZZ, CX, CY, CZ
View the matrix of each gate as

```
julia> cnot()            
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
```

Define a new gate. This is an alternative to defining the Matrix directly in the crcuit. It is helpful if a custom gate is used in multiple places or circuits

```
H(qubitOrder=missing) = gate([1 1; 1 -1]/sqrt(2), qubitOrder, "H")
```
The function paramter qubitOrder is optional while defining the gate but is nedded while using it in circuit. This is because the gate called without the qubitOrder returns the matrix representation of the gate.
There is a global variable "T" in the module that defines the precision. If you are planning to use higher precision, define gates as

```
H(qubitOrder=missing) = gate([T(1) T(1); T(1) T(-1)]/sqrt(T(2)), qubitOrder, "H")
```

## Importand Notes
This does not have much optimizations and all calculations are done in dense matrices.
Hence, this module is not preferable on high qubit count circuits.