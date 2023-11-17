module SimpleQuantumCircuit

using LinearAlgebra

export circuit, plotCircuit, gate, control, Ugate, Igate, X, Y,
Z,H,RX,RY,RZ,swap,cnot,XX,YY,ZZ,setType,
CX, CY, CZ

# Type definition that will be used for calculations across function in this module
T = Float64


"""
Convert string to Float64
"""
Float64(x::String) = parse(Float64,x)
#Float64(x::Number) = Base.Float64(x)


"""
Set global Type stored in variable T in the module.
This is used to define the precision of all computations by different functions of the module.
"""
function setType(type::Type)
    global T = type    
end

"""
Reverse of digits() function
    dig - digits array
    base - default base is 2
    returns integer
"""
function undigit(dig; base=2)
    (promoSum, promoBase) = promote(zero(eltype(dig)), base)
    multiplier = one(promoSum)
    for value in dig
        promoSum += value * multiplier
        multiplier *= promoBase
    end
    return Int(promoSum)
end

"""
Find index of matrix of gate based on index of gated Circuit
    n - total number of qubits
    oldIndex - gated Circuit Matrix index
    order - order of qubits that are used in gate
    returns index of matrix of gate
"""
function newIndex(n, oldIndex, order)
    digitizedIndex=digits(oldIndex, base=2, pad=n)
    digitizedNewIndex = digitizedIndex[order]
    undigit(digitizedNewIndex)
end

"""
expand the gate based on order of qubits and circuit
    n - total number of qubits
    gateMat - matrix representing gate
    qubitOrder - arrangement of circuit qubits that are used in gate
Example '''expandGate(3, CX, [1 3])'''
"""
function expandGate(n, gateMat, qubitOrder)
    circuitDim = 2^n
    # circuitMat = zeros(Complex{BigFloat}, circuitDim,circuitDim)
    circuitMat = zeros(Complex{T}, circuitDim,circuitDim)
    gateIndex = ones(Int, circuitDim)
    gateExcludedQubits = setdiff(1:n, qubitOrder)       # circuit qubits that are not part of the gate
    gateExcludedIndex = ones(Int, circuitDim)

    Threads.@threads for i = 1:circuitDim                               # calculates the indexes. This is done is a separate loop so that Indexes are calculated only once per matrix
      gateIndex[i] = newIndex(n, i-1, qubitOrder) + 1                   # new index in U for the gate matrix elements
      gateExcludedIndex[i] = newIndex(n, i-1, gateExcludedQubits) + 1   # value/spin of qubits not part of gate matrix elements
    end

    Threads.@threads for i = 1:circuitDim
      for j = 1:circuitDim
        if gateExcludedIndex[i] == gateExcludedIndex[j]                 # qubits that are not part of gate matrix sould have same value to apply gate elements.
          circuitMat[i,j] = gateMat[gateIndex[i],gateIndex[j]]
        end
      end
    end

    return circuitMat
end

"""
create a unitary evolution matrix from a quantum circuit
  n - number of qubits
  gateList - sequence of gates
Example circuit(5,(H,[1]),(CX,[1,2]))  
"""
function circuit(n, gateList...)
    plotCircuit(n, gateList...)
    circuitDim = 2^n
    U = Matrix{T}(I(circuitDim))
    for gate in gateList
        if T == BigFloat
            U = BigFloat_matrix_multiply(expandGate(n, gate[1], gate[2]),U)
        else
            U = expandGate(n, gate[1], gate[2])*U
        end
    end
    return U
end

"""
plot a circuit
    n - number of qubits, 
    gateList... - list of gates
"""
function plotCircuit(n, gateList...)
    StrVector = String[]

    # initialize
    for i in 1:n
        push!(StrVector, "   ")
        push!(StrVector, "$i--")
        push!(StrVector, "   ")
    end
    
    #prepare for each gate
    for gate in gateList
        for qubit in 1:n
          line = 3*(qubit-1)+1
          if qubit in gate[2]
            gateQubit = indexin(qubit, gate[2])[1]
            StrVector[line] = StrVector[line]*"┌──────┐"
            StrVector[line+1] = StrVector[line+1]*"┤$gateQubit "*rpad(gate[3], 4, " ")[1:4]*"├"
            StrVector[line+2] = StrVector[line+2]*"└──────┘"
          else
            StrVector[line] = StrVector[line]*"        "
            StrVector[line+1] = StrVector[line+1]*"--------"
            StrVector[line+2] = StrVector[line+2]*"        "
          end
        end
    end

    # print
    for i in 1:3*n
        println(StrVector[i])
    end   
end

"""
generic gate. If qubit Order is missing, it just return the gate matrix. 
"""
function gate(mat, qubitOrder=missing, label="UG")
    if qubitOrder === missing
        return mat
    else
        return (mat, qubitOrder, label)
    end
end

"""
control gate
  numControl - number of control qubits
  The qubits starting from first qubit are considered control qubits
  U - The matrix that is applied on rest of qubits when all the control qubits are 1.
"""
function control(numControl, U)
    nU::Int = log2(size(U)[1])
    n::Int = numControl + nU
    gateDim::Int = 2^n
    controlDim::Int = 2^numControl
    # gateMat = zeros(Complex{BigFloat}, gateDim,gateDim)
    gateMat = zeros(Complex{T}, gateDim,gateDim)

    for i = 1:gateDim
      for j = 1:gateDim
        if i == j              # If control qubits are all high
          gateMat[i,j] = 1.0
        end
        if i % controlDim == 0 && j % controlDim == 0
            x = Int(i / controlDim)
            y = Int(j / controlDim)
            gateMat[i,j] = U[x,y]
        end
      end
    end

    return gateMat
end

"""
Define own multithreaded matrix function to multiply BigFloat Complex Matrices
"""
function BigFloat_matrix_multiply(A::Matrix, B::Matrix)::Matrix
  rows_A, cols_A = size(A)
  rows_B, cols_B = size(B)

  # Ensure that the matrices can be multiplied
  if cols_A != rows_B
      error("Number of columns of A must be equal to number of rows of B for matrix multiplication.")
  end

  # Create a result matrix of zeros
  C = zeros(BigFloat, rows_A, cols_B)

  Threads.@threads for i = 1:rows_A
      for j = 1:cols_B
          for k = 1:cols_A
              C[i, j] += A[i, k] * B[k, j]
          end
      end
  end

  return C
end

"""
All defined gates
"""
Ugate(theta , phi, lambda, qubitOrder=missing) = begin theta = T(string(theta)); phi = T(string(phi)); lambda = T(string(lambda)); gate([cos(theta/T(2)) -exp(lambda*im)*sin(theta/T(2)); exp(phi*im)*sin(theta/T(2)) exp((phi+lambda)*im)*cos(theta/T(2))], qubitOrder, "U") end
Igate(qubitOrder=missing) = gate([T(1) T(0); T(0) T(1)], qubitOrder, "I")
X(qubitOrder=missing) = gate([T(0) T(1); T(1) T(0)], qubitOrder, "X")
Y(qubitOrder=missing) = gate([T(0) T(-1)*im; T(1)*im T(0)], qubitOrder, "Y")
Z(qubitOrder=missing) = gate([T(1) T(0); T(0) T(-1)], qubitOrder, "Z")
H(qubitOrder=missing) = gate([T(1) T(1); T(1) T(-1)]/sqrt(T(2)), qubitOrder, "H")
RX(theta, qubitOrder=missing) = begin theta = T(string(theta)); gate([cos(theta/T(2)) -im*sin(theta/T(2)); -im*sin(theta/T(2)) cos(theta/T(2))], qubitOrder, "RX")  end
RY(theta, qubitOrder=missing) = begin theta = T(string(theta)); gate([cos(theta/T(2)) -sin(theta/T(2)); sin(theta/T(2)) cos(theta/T(2))], qubitOrder, "RY") end
RZ(theta, qubitOrder=missing) = begin theta = T(string(theta)); gate([exp(-im*theta/T(2)) T(0); T(0) exp(im*theta/T(2))], qubitOrder, "RZ") end

swap(qubitOrder=missing) = gate([T(1) T(0) T(0) T(0); T(0) T(0) T(1) T(0); T(0) T(1) T(0) T(0); T(0) T(0) T(0) T(1)], qubitOrder, "SWAP")
cnot(qubitOrder=missing) = gate(control(1, X()), qubitOrder, "CNOT")
XX(theta, qubitOrder=missing) = begin theta = T(string(theta)); gate([cos(theta/T(2)) T(0) T(0) -im*sin(theta/T(2)); T(0) cos(theta/T(2))  -im*sin(theta/T(2)) T(0); T(0) -im*sin(theta/T(2)) cos(theta/T(2)) 0; -im*sin(theta/T(2)) T(0) T(0) cos(theta/T(2)) ], qubitOrder, "XX") end
YY(theta, qubitOrder=missing) = begin theta = T(string(theta)); gate([cos(theta/2) T(0) T(0) im*sin(theta/T(2)); 0 cos(theta/T(2))  -im*sin(theta/T(2)) T(0); T(0) -im*sin(theta/T(2)) cos(theta/T(2)) T(0); im*sin(theta/T(2)) T(0) T(0) cos(theta/T(2)) ], qubitOrder, "YY") end
ZZ(theta, qubitOrder=missing) = begin theta = T(string(theta)); gate([exp(-im*theta/T(2)) T(0) T(0) T(0); T(0) exp(im*theta/T(2))  T(0) T(0); T(0) T(0) exp(im*theta/T(2)) T(0); T(0) T(0) T(0) exp(-im*theta/T(2)) ], qubitOrder, "ZZ") end

CX(qubitOrder=missing) = gate(control(1,X()), qubitOrder, "CX")
CY(qubitOrder=missing) = gate(control(1,Y()), qubitOrder, "CY")
CZ(qubitOrder=missing) = gate(control(1,Z()), qubitOrder, "CZ")

end
