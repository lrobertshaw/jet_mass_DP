samples=("LH_20_qq" "LH_30_qq" "LH_40_qq" "LH_50_qq" "LH_60_qq" "LH_70_qq" "LH_80_qq" "LH_90_qq" "MX_350_LH_50_qq")
# samples=("MinBias" "TT_boosted" "LH_20_qq" "LH_30_qq" "LH_40_qq" "LH_50_qq" "LH_60_qq" "LH_70_qq" "LH_80_qq" "LH_90_qq" "MX_350_LH_50_qq")
base_cfg="configs/V49nano_AR25/rate_table/step1and2_cfg.yml"

for sample in "${samples[@]}"; do
    cfg="configs/V49nano_AR25/rate_table/${sample}_cfg.yml"
    sed -e "s/^sample:.*/sample: $sample/" \
        -e "s/^table_fname:.*/table_fname: ${sample}_RateTable/" \
        "$base_cfg" > "$cfg"

    rate_table "$cfg" &  # run in background
done

wait  # wait for all jobs to finish