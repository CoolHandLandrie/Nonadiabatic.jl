using Distributed
using ClusterManagers
using DelimitedFiles

function __init__()
    include("/lustre/isaac/proj/UTK0015/dak/NonAdiabatic/src/FSSH.jl")
end


# Function to create Gaussian16 input file
function create_g16_input(molecule::String, method::String, basis_set::String, excited_states::Int)
    input_content = """
    %chk=$(molecule).chk
    %mem=145GB
    %nprocshared=40
    #p $method/$basis_set empiricaldispersion=gd3bj SCF=(XQC) guess=Read Geom=AllCheck force TD=(Nstates=$excited_states,NAC)

    $molecule excited state calculation

    """
    open("$(molecule).com", "w") do file
        write(file, input_content)
    end
    return "$(molecule).com"
end

# Main workflow
function prepare_excited_state_calculation(molecule::String, method::String, basis_set::String, excited_states::Int)
    input_file = create_g16_input(molecule, method, basis_set, excited_states)
    println("Input file created: $input_file")
    println("You can now submit this job manually using your preferred method.")
    return input_file
end

# Function to assign masses based on atom type
function assign_mass(atom_type::String)
    masses = Dict(
        "H" => 1.008,
        "C" => 12.011,
        "N" => 14.007,
        "S" => 32.065,
        "Se" => 78.971,
        "Te" => 127.60,
        "Fe" => 55.845,
        "Cu" => 63.546,
        "Ag" => 107.868,
        "Au" => 196.967,
        # Add more elements as required
    )
    return get(masses, atom_type, NaN)  # Returns NaN if atom type not found
end

function extract_energies(gaussian_output::String)
    energies = Float64[]
    open(gaussian_output, "r") do file
        for line in eachline(file)
            if occursin("Excited State", line)
                energy = parse(Float64, split(line)[4])
                push!(energies, energy)
            end
        end
    end
    return [0.0; energies]  # Add ground state energy (0.0) at the beginning
end

function extract_couplings(gaussian_output::String)
    couplings = Dict{Tuple{Int,Int}, Float64}()
    current_states = (0, 0)
    open(gaussian_output, "r") do file
        for line in eachline(file)
            if occursin("NACMEs between states", line)
                current_states = Tuple(parse.(Int, split(line)[end-1:end]))
            elseif occursin("dE/dX", line) && current_states != (0, 0)
                coupling = parse(Float64, split(line)[end])
                couplings[current_states] = coupling
            end
        end
    end
    return couplings
end

function extract_forces(gaussian_output::String)
    forces = Vector{Float64}[]
    reading_forces = false
    open(gaussian_output, "r") do file
        for line in eachline(file)
            if occursin("Forces (Hartrees/Bohr)", line)
                reading_forces = true
                continue
            end
            if reading_forces
                if occursin("-------------------", line)
                    reading_forces = false
                    break
                end
                force = parse.(Float64, split(line)[3:5])
                push!(forces, force)
            end
        end
    end
    return forces
end

function julia_main()::Cint
    # Run the calculation
    molecule = "your_complex_molecule"
    method = "ωB97X-D"
    basis_set = "def2-TZVPP"
    excited_states = 10

    # Prepare Gaussian16 input
    input_file = prepare_excited_state_calculation(molecule, method, basis_set, excited_states)

    # After running Gaussian16 and extracting necessary data, set up FSSH simulation

    # Define atom types in your molecule (example, adjust as needed)
    atom_types = ["Au", "Ag", "Cu", "N", "C", "C", "C", "H", "H", "H", "S", "Se", "Te"]
    masses = [assign_mass(atom) for atom in atom_types]

    # Define parameters for FSSH
    initial_state = 1  # Start in ground state
    num_states = excited_states
    num_steps = 1000
    dt = 0.1  # Time step in fs

    # Run FSSH simulation
    gaussian_output = "path/to/your/gaussian_output_file.log"  # Replace with actual path
    trajectory = run_fssh(
        initial_state,
        num_states,
        num_steps,
        dt,
        extract_energies,
        extract_couplings,
        extract_forces,
        masses,
        gaussian_output
    )

    # Add any output or result processing here
    println("FSSH simulation completed. Trajectory data:")
    println(trajectory)

    return 0  # Return 0 to indicate successful execution
end

# This line is necessary for creating a standalone executable
Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    return julia_main()
end

# If you want to be able to run this script directly in Julia, uncomment the following line:
# if abspath(PROGRAM_FILE) == @__FILE__
#     julia_main()
# end

