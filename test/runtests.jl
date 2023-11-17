using SimpleQuantumCircuit
using Test

@testset "SimpleQuantumCircuit.jl" begin
    @test SimpleQuantumCircuit.X() != 0       #Checks Module is compilable
    @test SimpleQuantumCircuit.X() == [0 1; 1 0]  # Check X gate matrix
    @test SimpleQuantumCircuit.X([1])[1] == [0 1; 1 0]    # check gate() function
    @test SimpleQuantumCircuit.X([1])[2] == [1]    # check gate() function
    @test SimpleQuantumCircuit.X([1])[3] == "X"    # check gate() function
    @test SimpleQuantumCircuit.cnot() == [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]  # Check control gate
    @test SimpleQuantumCircuit.circuit(2, cnot([1,2]) ) == [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]  # Check building circuits

    # check higher precision
    setType(BigFloat)
    setprecision(128)
    @test typeof(H()) == Matrix{BigFloat}

    # check higher percision circuits
    @test SimpleQuantumCircuit.circuit(2, cnot([1,2]) ) == [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0] 
end
