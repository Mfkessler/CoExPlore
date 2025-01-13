#!/usr/bin/env python3
"""
simulate_testdata.py

Generates simulated RNA-seq test data:
1) A count matrix per species (pseudo-TMM values)
2) A sample mapping (one per species)
3) GO-IDs and InterPro-IDs (IPR) per transcript
4) An orthogroup file (N0.tsv) for n>1 species with Orthofinder-like structure
"""

import numpy as np
import pandas as pd
import argparse
import os
import random
from collections import defaultdict


def generate_sample_mapping(species_id, traits, samples_per_trait):
    """
    Generates a pd.DataFrame with columns:
      - sample
      - Testgroup (corresponds to "Trait" in this case)
    e.g.
    sample    Testgroup
    S1_A_1    Trait A
    S1_A_2    Trait A
    ...
    """
    rows = []
    for trait in traits:
        for rep in range(samples_per_trait):
            sample_name = f"{species_id}_{trait}_{rep+1}"
            rows.append((sample_name, trait))
    df = pd.DataFrame(rows, columns=["sample", "Testgroup"])
    return df


def generate_modules(num_genes, num_modules, module_size_mean, module_size_sd):
    """
    Randomly divides 'num_genes' into 'num_modules' + rest module.
    Returns an array of length num_genes, where each gene has the module ID (-1 for "no module"/noise).
    """
    # Normally distributed module sizes (round positively, then adjust)
    # Alternatively, you can use random.exponential etc.
    mod_sizes = np.abs(np.round(np.random.normal(module_size_mean, module_size_sd, num_modules))).astype(int)
    # If too few or too many genes are covered, normalize a bit
    total_mod = np.sum(mod_sizes)
    if total_mod > num_genes:
        # If the sum is greater than num_genes, shorten some modules
        scale_factor = num_genes / total_mod
        mod_sizes = np.floor(mod_sizes * scale_factor).astype(int)
        total_mod = np.sum(mod_sizes)
    # Random order
    gene_indices = np.arange(num_genes)
    np.random.shuffle(gene_indices)
    module_assignment = np.full(num_genes, -1, dtype=int)  # -1 = no module
    start_idx = 0
    for i, size in enumerate(mod_sizes):
        end_idx = start_idx + size
        if end_idx > num_genes:
            end_idx = num_genes
        module_assignment[gene_indices[start_idx:end_idx]] = i
        start_idx = end_idx
    return module_assignment


def generate_counts(num_genes, samples, module_assignment, cor_strength=0.8, noise_level=0.2):
    """
    Generates a "pseudo-TMM" matrix with shape=(num_genes, len(samples)).

    - For each module, create a "base profile" for each sample (e.g. normally distributed),
      which is strongly correlated among the genes of the same module.
    - noise_level controls the amount of noise
    - cor_strength controls how similar the genes in a module are
    """
    num_samples = len(samples)
    # Initialize count matrix
    data = np.zeros((num_genes, num_samples), dtype=float)

    # Create a "profile" per module
    unique_modules = np.unique(module_assignment[module_assignment >= 0])  # without -1
    module_profiles = dict()

    for mod_id in unique_modules:
        # For each sample, a module base value
        # e.g. draw normal(5,2) => "average expression"
        base_profile = np.random.normal(loc=5.0, scale=2.0, size=num_samples)
        module_profiles[mod_id] = base_profile
    
    # Then for each gene, if it is in module X, generate
    # Expression ~ base_profile + correlated_component + noise
    for g in range(num_genes):
        m = module_assignment[g]
        if m == -1:
            # Noise-only
            gene_profile = np.random.normal(loc=2.0, scale=3.0, size=num_samples)
        else:
            base = module_profiles[m]
            # correlated component
            correlated_part = np.random.normal(0, (1 - cor_strength), size=num_samples)
            # Sum and noise
            gene_profile = base + cor_strength*correlated_part + noise_level*np.random.normal(0, 1, size=num_samples)
        # exponentiate to get TMM-like values
        # (or use the value directly, depending on preference)
        gene_profile = np.exp(gene_profile)  # exponential => can vary widely
        data[g, :] = gene_profile
    
    # final: could transform the data to log2 again
    # but leave it as "exponential" => pseudo-TMM
    # Optional: Scale the data per sample to introduce "library size" variation
    lib_sizes = np.random.uniform(0.8, 1.2, size=num_samples)
    data = data * lib_sizes  # slight scaling per column

    return data


def generate_go_ipr(num_genes, go_per_gene=3, ipr_per_gene=2,
                    go_range=(0, 9999999), ipr_range=(0, 999999),
                    ipr_desc_pool=None):
    """
    Generates a list of GO-IDs and InterPro-IDs + Description per gene.
    ipr_desc_pool: List of possible descriptions (e.g. ["Desc1", "Desc2", ...])
    """
    # Create ipr_desc_pool generically if needed
    if ipr_desc_pool is None:
        ipr_desc_pool = [f"Desc{i}" for i in range(1, 2001)]
    
    go_data = []
    ipr_data = []
    for g in range(num_genes):
        # GO
        # Randomly select go_per_gene GOs from the range
        gos = []
        for _ in range(go_per_gene):
            val = np.random.randint(go_range[0], go_range[1] + 1)
            gos.append(f"GO:{val:07d}")  # => GO:0001234
        # IPR
        iprs = []
        for _ in range(ipr_per_gene):
            val = np.random.randint(ipr_range[0], ipr_range[1] + 1)
            ipr_id = f"IPR{val:06d}"
            desc = random.choice(ipr_desc_pool)
            iprs.append((ipr_id, desc))
        go_data.append(gos)
        ipr_data.append(iprs)
    return go_data, ipr_data


def write_count_matrix(file_path, data, samples, transcript_ids):
    """
    Format:
    Header: [tab] sample1 sample2 ...
    Rows: transcript_id <tab> val1 val2 ...
    """
    # Create DataFrame
    df = pd.DataFrame(data, index=transcript_ids, columns=samples)
    # Index as first column
    df.index.name = ""
    # Output DataFrame via to_csv (tab-sep)
    df.to_csv(file_path, sep="\t", index=True, header=True)


def write_sample_mapping(file_path, sample_mapping_df):
    """
    sample_mapping_df with columns: sample, Testgroup
    """
    sample_mapping_df.to_csv(file_path, sep="\t", index=False)


def write_go_file(file_path, transcript_ids, go_data):
    """
    Format:
    Transcript1  GO:0000000 GO:0001234
    ...
    -> transcript_id \t GO_1 GO_2 ...
    """
    with open(file_path, "w") as f:
        for t_id, gos in zip(transcript_ids, go_data):
            gos_str = " ".join(gos)
            f.write(f"{t_id}\t{gos_str}\n")


def write_ipr_file(file_path, transcript_ids, ipr_data):
    """
    Format:
    Transcript1  IPR003593 desc
    Transcript1  IPR003959 desc
    ...
    => transcript_id \t ipr_id \t description
    """
    with open(file_path, "w") as f:
        for t_id, iprs in zip(transcript_ids, ipr_data):
            for (ipr_id, desc) in iprs:
                f.write(f"{t_id}\t{ipr_id}\t{desc}\n")


def generate_orthogroups(num_species, transcript_ids_per_species, coverage=0.7):
    """
    Generates pseudo-orthogroups. coverage=0.7 => 70% of transcripts have orthologs in other species.
    Generate HOG IDs ("HOG1", "HOG2", ...), OG IDs ("OG1", "OG2", ...).
    Format: "HOG     OG      Gene Tree Parent Clade  S1 ... S2 ... S3 ..."
    """
    # All transcripts in a pool
    # Form random orthogroups => e.g. for each HOG 1..n pick 1-3 genes per species (depending on coverage)
    # Create e.g. 1 orthogroup per 10 transcripts (or so).
    # The input coverage controls: How many transcripts should be in orthogroups?
    og_rows = []
    hog_id = 1
    og_id = 1
    # Flatten
    total_tids = []
    for sp_idx in range(num_species):
        for t_id in transcript_ids_per_species[sp_idx]:
            total_tids.append((sp_idx, t_id))
    
    random.shuffle(total_tids)
    # put coverage * len(total_tids) transcripts into orthogroups
    keep_count = int(coverage * len(total_tids))
    used = set()  # Set of (sp_idx, t_id) that have already been assigned
    i = 0
    while i < keep_count:
        # Create a new orthogroup
        # Choose 0-1 transcript per species? or 1-n per species?
        # "0 or 1 per species" (simplified)
        group_members = defaultdict(list)
        for sp_idx in range(num_species):
            # 30% chance to take a transcript from this species
            if random.random() < 0.6:
                # take the next free entry
                found = False
                while i < len(total_tids):
                    if total_tids[i] not in used and total_tids[i][0] == sp_idx:
                        group_members[sp_idx].append(total_tids[i][1])
                        used.add(total_tids[i])
                        i += 1
                        found = True
                        break
                    i += 1
                if not found:
                    # Species has no free transcripts left
                    pass
        
        # Build row
        # HOG, OG, Gene Tree Parent Clade => doesn't matter
        row_dict = {
            "HOG": f"HOG{hog_id}",
            "OG": f"OG{og_id}",
            "Gene Tree Parent Clade": f"CladeX{hog_id}"
        }
        # Now S1..S{num_species} columns
        for sp_idx in range(num_species):
            sp_name = f"S{sp_idx+1}"
            members_str = ",".join(group_members[sp_idx]) if sp_idx in group_members else ""
            row_dict[sp_name] = members_str
        og_rows.append(row_dict)
        hog_id += 1
        og_id += 1
    
    # Now build a DataFrame
    columns = ["HOG", "OG", "Gene Tree Parent Clade"]
    for sp_idx in range(num_species):
        columns.append(f"S{sp_idx+1}")
    df = pd.DataFrame(og_rows, columns=columns)
    return df


def write_n0_file(file_path, orthogroup_df):
    """
    Format as desired:
    HOG     OG      Gene Tree Parent Clade  S1      S2      ... SN
    """
    orthogroup_df.to_csv(file_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description="Simulate bioinformatics test data for any number of species.")
    parser.add_argument("--output_dir", type=str, default="../input/testdata_out",
                        help="Target directory for output files.")
    parser.add_argument("--num_species", type=int, default=2,
                        help="Number of species (S1..S{n}).")
    parser.add_argument("--num_transcripts", type=int, default=5000,
                        help="Number of transcripts per species.")
    parser.add_argument("--num_traits", type=int, default=22,
                        help="Number of traits (max. 22).")
    parser.add_argument("--samples_per_trait", type=int, default=5,
                        help="Number of samples per trait.")
    parser.add_argument("--num_modules", type=int, default=10,
                        help="Number of correlated modules.")
    parser.add_argument("--module_size_mean", type=float, default=300.0,
                        help="Average size of a module (normally distributed).")
    parser.add_argument("--module_size_sd", type=float, default=100.0,
                        help="Standard deviation of module sizes.")
    parser.add_argument("--cor_strength", type=float, default=0.8,
                        help="Correlation within a module (0..1).")
    parser.add_argument("--noise_level", type=float, default=0.2,
                        help="Noise level.")
    parser.add_argument("--go_per_gene", type=int, default=2,
                        help="How many GO-IDs per gene/transcript.")
    parser.add_argument("--ipr_per_gene", type=int, default=2,
                        help="How many IPR-IDs per gene/transcript.")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility.")
    args = parser.parse_args()

    # Seed
    np.random.seed(args.seed)
    random.seed(args.seed)

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Trait names (up to 22 tissues possible)
    possible_traits = [
        "Bud stage 1",
        "Bud stage 2",
        "Bud stage 3",
        "Bud stage 4",
        "Sepals at anthesis",
        "Petals at anthesis",
        "Stamens at anthesis",
        "Gynoecia at anthesis",
        "Shoot apex",
        "Petal early stage",
        "Petal mid stage",
        "Petal late stage",
        "Young fruits",
        "Mid-stage fruit",
        "Seeds 5 dap",
        "Root",
        "Young leaf",
        "Mature leaf",
        "Seedling",
        "Mature petal nectary part",
        "Mature petal no nectary part",
        "Non-spurred Petal"
    ]

    # Plausibility check
    if args.num_traits > len(possible_traits):
        raise ValueError("num_traits must be max. 22.")

    # Trim the list to the desired number of traits
    traits = possible_traits[:args.num_traits]


    # For each species generate:
    # - transcript IDs
    # - sample mapping
    # - modules
    # - counts
    # - GO, IPR

    all_transcript_ids = []
    for sp_idx in range(args.num_species):
        species_name = f"S{sp_idx+1}"  # e.g. S1, S2 etc.
        print(f"Creating species {species_name}...")

        # Transcript IDs
        transcript_ids = [f"Transcript{i+1}" for i in range(args.num_transcripts)]
        all_transcript_ids.append(transcript_ids)

        # Sample mapping
        sample_mapping_df = generate_sample_mapping(species_name, traits, args.samples_per_trait)
        # Save
        sample_mapping_path = os.path.join(args.output_dir, f"{species_name}.sample_mapping.tsv")
        write_sample_mapping(sample_mapping_path, sample_mapping_df)

        # Module assignment
        module_assignment = generate_modules(
            args.num_transcripts,
            args.num_modules,
            args.module_size_mean,
            args.module_size_sd
        )

        # Counts
        samples = sample_mapping_df["sample"].tolist()
        count_data = generate_counts(
            args.num_transcripts,
            samples,
            module_assignment,
            cor_strength=args.cor_strength,
            noise_level=args.noise_level
        )
        # Write
        count_matrix_path = os.path.join(args.output_dir, f"{species_name}.isoform.TMM.matrix")
        write_count_matrix(count_matrix_path, count_data, samples, transcript_ids)

        # GO, IPR
        go_data, ipr_data = generate_go_ipr(
            args.num_transcripts,
            go_per_gene=args.go_per_gene,
            ipr_per_gene=args.ipr_per_gene
        )
        # Write
        go_file_path = os.path.join(args.output_dir, f"{species_name}.goid.tsv")
        write_go_file(go_file_path, transcript_ids, go_data)
        ipr_file_path = os.path.join(args.output_dir, f"{species_name}.iprid.tsv")
        write_ipr_file(ipr_file_path, transcript_ids, ipr_data)

    # Orthogroups only if num_species > 1
    if args.num_species > 1:
        # Generate Orthogroup-N0.tsv
        print("Creating orthogroups (N0.tsv) ...")
        og_df = generate_orthogroups(args.num_species, all_transcript_ids, coverage=0.7)
        n0_path = os.path.join(args.output_dir, "N0.tsv")
        write_n0_file(n0_path, og_df)

    print("Done! All files are in", args.output_dir)


if __name__ == "__main__":
    main()
