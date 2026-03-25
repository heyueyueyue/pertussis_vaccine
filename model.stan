functions {
    // This function calculates the changes in compartments for the SEIRS model dynamically for each age group
    vector disease_dynamics(
        int N,                     // Number of age groups
        real B,                    // Birth rate
        real beta,                 // Transmission rate
        real sigma,                // Rate from exposed to infectious
        real gamma,                // Recovery rate
        real v,                    // Vaccination rate
        real lambda,               // Vaccine waning rate
        real alpha,                // Loss of immunity rate
        vector mu,                 // Aging rate between groups
        vector delta,              // Death rate
        vector S,                  // Susceptibles
        vector V,                  // Vaccinated
        vector E,                  // Exposed
        vector I,                  // Infectious
        vector R,                  // Recovered
        matrix contact_matrix      // Contact matrix between groups
    ) 
    {
        vector[N] dS = rep_vector(0, N);
        vector[N] dV = rep_vector(0, N);
        vector[N] dE = rep_vector(0, N);
        vector[N] dI = rep_vector(0, N);
        vector[N] dR = rep_vector(0, N);
        
        vector[N] lambda_infection = beta * (contact_matrix * I);
        vector[N] sigma_rate = rep_vector(1 / sigma, N);
        vector[N] gamma_rate = rep_vector(1 / gamma, N);
        vector[N] alpha_rate = rep_vector(1 / alpha, N);
        
        // Vaccine dynamics
        dV -= V .* lambda + delta .* V;
        dV[3:5] += v * S[3:5]; // Vaccinating susceptible individuals in specific age groups
        dV[8] += v * S[8];
        
        // Dynamics of susceptibles, exposed, infectious, and recovered
        dS += alpha_rate .* R - lambda_infection .* S - delta .* S;
        dS[1] += B;  // Adding newborns to the susceptible pool
        dS[3:5] -= v * S[3:5];
        dS[8] -= v * S[8];
        
        dE += lambda_infection .* S - sigma_rate .* E - delta .* E;
        dI += sigma_rate .* E - gamma_rate .* I - delta .* I;
        dR += gamma_rate .* I - alpha_rate .* R - delta .* R;
        
        // Age group transitions
        for (age in 1:(N-1)) {
            dS[age + 1] += mu[age] * S[age];
            dV[age + 1] += mu[age] * V[age];
            dE[age + 1] += mu[age] * E[age];
            dI[age + 1] += mu[age] * I[age];
            dR[age + 1] += mu[age] * R[age];

            dS[age] -= mu[age] * S[age];
            dV[age] -= mu[age] * V[age];
            dE[age] -= mu[age] * E[age];
            dI[age] -= mu[age] * I[age];
            dR[age] -= mu[age] * R[age];
        }

        return to_vector(append_row(append_row(append_row(append_row(dS, dV), dE), dI), dR));
    }
}

data {
    int<lower=1> N;  // Number of age groups
    matrix[N, N] contact_matrix;  // Contact matrix
    real<lower=0> v;  // Vaccination rate
    real<lower=0> lambda;  // Vaccine waning rate
    real<lower=0> B;  // Birth rate
    real<lower=0> N_pop;  // Total population
    real<lower=0> beta;  // Transmission rate
    vector<lower=0>[N] mu;  // Aging rate between groups
    vector<lower=0>[N] delta;  // Death rate
    int<lower=1> T;  // Number of days observed
    int<lower=0> observed_new_cases[N, T];  // New cases per day per age group
}

parameters {
    real<lower=0> sigma;  // Rate from exposed to infectious
    real<lower=0> gamma;  // Recovery rate
    real<lower=0> alpha;  // Loss of immunity rate
    vector<lower=0>[N] S;  // Initial number of susceptibles
    vector<lower=0>[N] V;  // Initial number of vaccinated
    vector<lower=0>[N] E;  // Initial number of exposed
    vector<lower=0>[N] I;  // Initial number of infectious
    vector<lower=0>[N] R;  // Initial number of recovered
    matrix[T, 5 * N] y;  // All model states over time
}

model {
    sigma ~ gamma(2, 0.2);  // Prior for the rate from exposed to infectious
    gamma ~ gamma(2, 0.1);  // Prior for the recovery rate
    alpha ~ lognormal(log(8 * 365), 0.2);  // Prior for the loss of immunity rate

    S + V + E + I + R ~ dirichlet(rep_vector(1, 5 * N));  // Prior for initial state proportions

    vector[5 * N] dydt;
    dydt = disease_dynamics(N, B, beta, sigma, gamma, v, lambda, alpha, mu, delta, S, V, E, I, R, contact_matrix);
    y[1] ~ normal(dydt, 0.1);  // Dynamics for the first day
    observed_new_cases[:, 1] ~ poisson(y[1][3*N+1:4*N]);  // Observing new cases for day one

     for (t in 2:T) {
         // Convert matrix rows to vectors
         vector[5 * N] y_prev = to_vector(y[t-1]);
     
         // Compute dynamics using the converted vectors
         dydt = disease_dynamics(N, B, beta, sigma, gamma, v, lambda, alpha, mu, delta,
                                  to_vector(y_prev[1:N]), to_vector(y_prev[N+1:2*N]),
                                  to_vector(y_prev[2*N+1:3*N]), to_vector(y_prev[3*N+1:4*N]),
                                  to_vector(y_prev[4*N+1:5*N]), contact_matrix);
     
         // Ensure that both y_prev and dydt are vectors for the addition to be valid
         y[t] ~ normal(y_prev + dydt, 0.1);  // Dynamics for subsequent days
     
         // Convert the relevant segment of y[t] back to vector if necessary for correct dimension matching
         observed_new_cases[:, t] ~ poisson(to_vector(y[t][(3*N+1):(4*N)]));  // Observing new cases per day
     }
}
