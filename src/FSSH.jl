# Fewest Switches Surface Hopping function
function fewest_switches_surface_hopping(
    initial_state::Int,
    num_states::Int,
    num_particles::Int,
    time_steps::Int,
    dt::Float64,
    energies::Function,
    couplings::Function,
    forces::Function,
    masses::Vector{Float64},
    initial_positions::Vector{Vector{Float64}},
    initial_velocities::Vector{Vector{Float64}}
)
    # Initialize variables
    current_state = initial_state
    amplitudes = zeros(ComplexF64, num_states)
    amplitudes[initial_state] = 1.0
    
    trajectory = []
    r = initial_positions
    v = initial_velocities
    
    # Quadrature rule
    rule = GaussLegendre(5)  # 5-point Gauss-Legendre quadrature
    
    for t in 1:time_steps
        # Get energies, couplings, and forces at current time step
        E = energies(t * dt)
        V = couplings(t * dt)
        F = forces(t * dt, current_state)
        
        # Propagate amplitudes
        H = E + V
        U = exp(-im * H * dt)
        amplitudes = U * amplitudes
        
        # Calculate hopping probabilities
        probs = zeros(num_states)
        for j in 1:num_states
            if j != current_state
                probs[j] = 2 * real(conj(amplitudes[current_state]) * amplitudes[j] * V[current_state, j]) * dt / abs(amplitudes[current_state])^2
            end
        end
        
        # Attempt hop
        if rand() < sum(probs)
            new_state = sample(1:num_states, Weights(probs))
            if new_state != current_state
                # Perform hop (adjust velocities, etc.)
                current_state = new_state
            end
        end
        
        # Update nuclear positions and velocities using quadrature
        for i in 1:num_particles
            # Update position
            r[i] += quadrature(rule, s -> v[i], t*dt, (t+1)*dt)
            
            # Update velocity
            v[i] += quadrature(rule, s -> F[i] / masses[i], t*dt, (t+1)*dt)
        end
        
        push!(trajectory, (t * dt, current_state, copy(amplitudes), copy(r), copy(v)))
    end
    
    return trajectory
end