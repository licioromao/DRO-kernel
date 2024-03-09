function noise = generate_noise_LTI(mean, chol_cov)
    
    noise = mean + chol_cov*randn(2,1);
    
end