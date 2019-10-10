#!/ bin / bash
#PBS - N Unfolding
#PBS - o output.out
#PBS - e output.err
#module load gcc / 4.8.5

module load gcc / 6.3.0 / ion / home / hect / Preprocessing / PreProcLenny / bin / pCT_Preprocessing / ion / pCT_data /
    raw_data / CPC_2016_08_13 / HIT_Dec / LinePair_0034_Cont_000.dat &wait
