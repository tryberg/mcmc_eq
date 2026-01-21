#!/bin/bash

# this is a readable version of the awk script within mkSynthetics courtesy of Claude

# Add gaussian noise to seismic arrival times
# Noise levels by quality code:
#   P-wave: P0=0.05s, P1=0.10s, P2=0.15s, P3=0.20s
#   S-wave: S0=0.175s, S1=0.225s, S2=0.275s, S3=0.325s

# Configuration
RMS_NOISE=0.03        # Base RMS noise level
RANDOM_SEED=33        # Random seed for reproducibility
RANDOM_SEED=$RANDOM   
INPUT_FILE="synths_wo_noise"
OUTPUT_FILE="synths_with_noise"

# Generate noisy data using Box-Muller transform for Gaussian noise
awk -v rms="$RMS_NOISE" -v seed="$RANDOM_SEED" '
BEGIN {
    # Initialize random seed
    srand(seed)
}
{
    # Generate Gaussian random number using Box-Muller transform
    do {
        x1 = 2.0 * rand() - 1    # Random number in [-1, 1]
        x2 = 2.0 * rand() - 1    # Random number in [-1, 1]
        s = x1*x1 + x2*x2        # Sum of squares
    } while ((s > 1) || (s == 0))  # Reject if outside unit circle or zero
    
    # Box-Muller transformation to get Gaussian distribution
    gaussian_noise = x1 * sqrt(-2.0 * log(s) / s)
    
    # Check if line is a comment/header
    if ($1 == "#") {
        print $0
    } else {
        # Determine if this is an S-wave (affects noise scaling)
        is_s_wave = 0
        if ($3 == "S") {
            is_s_wave = 1
        }
        
        # Calculate noise scaling factor based on quality code ($8)
        # Quality codes: 0, 1, 2, 3 (lower is better)
        # S-waves get 2.5x more noise than P-waves
        noise_scale = ((($8 + 1) + 2.5 * is_s_wave) / 4) * 2
        
        # Add scaled Gaussian noise to arrival time ($7)
        noisy_time = $7 + (gaussian_noise * rms * noise_scale)
        
        # Output: station phase_type arrival_time quality_code ...
        print $1, $2, $3, $4, $5, $6, noisy_time, $8
    }
}' "$INPUT_FILE" > "$OUTPUT_FILE"

# Calculate and display noise levels
p0=$(echo "$RMS_NOISE * ((0+1+2.5*0)/4) * 2" | bc -l)
p1=$(echo "$RMS_NOISE * ((1+1+2.5*0)/4) * 2" | bc -l)
p2=$(echo "$RMS_NOISE * ((2+1+2.5*0)/4) * 2" | bc -l)
p3=$(echo "$RMS_NOISE * ((3+1+2.5*0)/4) * 2" | bc -l)
s0=$(echo "$RMS_NOISE * ((0+1+2.5*1)/4) * 2" | bc -l)
s1=$(echo "$RMS_NOISE * ((1+1+2.5*1)/4) * 2" | bc -l)
s2=$(echo "$RMS_NOISE * ((2+1+2.5*1)/4) * 2" | bc -l)
s3=$(echo "$RMS_NOISE * ((3+1+2.5*1)/4) * 2" | bc -l)

echo "Added Gaussian noise to $INPUT_FILE -> $OUTPUT_FILE"
echo "RMS noise level: $RMS_NOISE"
echo "Random seed: $RANDOM_SEED"
echo ""
echo "Noise levels applied:"
printf "  P0: %.4f s\n" $p0
printf "  P1: %.4f s\n" $p1
printf "  P2: %.4f s\n" $p2
printf "  P3: %.4f s\n" $p3
printf "  S0: %.4f s\n" $s0
printf "  S1: %.4f s\n" $s1
printf "  S2: %.4f s\n" $s2
printf "  S3: %.4f s\n" $s3
