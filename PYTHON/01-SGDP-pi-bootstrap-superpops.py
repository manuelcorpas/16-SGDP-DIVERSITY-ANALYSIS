import pandas as pd
import allel
import numpy as np
import matplotlib.pyplot as plt
import sys
import re

# --- Paths ---
vcf_path = 'DATA/SGDP/cteam_extended.v4.PS2_phase.public.chr22.vcf.gz'
metadata_path = 'DATA/SGDP/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt'

# --- Modern global population representation by ancestry categories (2022 UN estimates) ---
# Data from GWAS Diversity Monitor - represents current global population distribution
modern_population_representation = {
    'European': 9.4,                              # Europe + European ancestry populations
    'Asian': 45.1,                                # East Asian (20.8%) + South Asian (24.3%)
    'African': 17.2,                              # Sub-Saharan Africa
    'African American or Afro-Caribbean': 1.0,    # African diaspora populations
    'Hispanic or Latin American': 8.1,            # Latin America and Hispanic populations
    'Other/Mixed': 18.2                           # Rest of world populations
}

# --- Load SGDP metadata ---
print("Loading SGDP metadata...")

# Try different encodings to handle potential encoding issues
encodings_to_try = ['utf-8', 'latin-1', 'ISO-8859-1', 'cp1252', 'utf-16']

metadata = None
successful_encoding = None

for encoding in encodings_to_try:
    try:
        print(f"Trying encoding: {encoding}")
        # First, let's check the file structure
        with open(metadata_path, 'r', encoding=encoding) as f:
            header_line = f.readline().strip()
            print(f"Header line: {header_line}")
        
        # Try tab-separated first (based on the header format)
        metadata = pd.read_csv(metadata_path, sep="\t", encoding=encoding)
        print(f"Successfully loaded metadata with {encoding} encoding")
        print(f"Shape: {metadata.shape}")
        print(f"Columns: {metadata.columns.tolist()}")
        successful_encoding = encoding
        break
        
    except Exception as e:
        print(f"Error with {encoding}: {e}")
        continue

if metadata is None:
    # If all encodings fail, try with error handling
    try:
        metadata = pd.read_csv(metadata_path, sep="\t", encoding='utf-8', errors='ignore')
        print("Loaded with UTF-8 ignoring errors")
    except Exception as e:
        print(f"Final attempt failed: {e}")
        raise RuntimeError("Could not load metadata file with any encoding")

print(f"First few rows:\n{metadata.head()}")

# Let's examine what populations are actually in the data
print("\n" + "="*60)
print("EXAMINING ACTUAL POPULATIONS IN SGDP METADATA")
print("="*60)

# Clean up column names
metadata.columns = [col.lstrip('#') for col in metadata.columns]

# Find the population column
pop_col = 'Population_ID'
if pop_col in metadata.columns:
    unique_pops = sorted(metadata[pop_col].unique())
    print(f"Total unique populations: {len(unique_pops)}")
    print("All populations in SGDP:")
    for i, pop in enumerate(unique_pops):
        print(f"{i+1:3d}. {pop}")
else:
    print("Could not find Population_ID column")
    print(f"Available columns: {metadata.columns.tolist()}")

# --- Map SGDP populations to ancestry categories ---
# Create mapping based on actual SGDP population names

sgdp_to_ancestry_category = {}

# African populations (Sub-Saharan Africa and North Africa)
african_patterns = [
    # Sub-Saharan Africa
    'San', 'Ju_hoan', 'Yoruba', 'Mbuti', 'Biaka', 'Aka', 'Baka', 
    'Mandenka', 'Dinka', 'Hadza', 'Sandawe', 'Luhya', 'Bantu', 
    'BantuSouthAfrica', 'BantuKenya', 'BantuTswana', 'BantuHerero',
    'Kongo', 'Esan', 'Lemande', 'Laka', 'Khomani_San', 'Masai',
    'Somali', 'Igbo', 'Kikuyu', 'Fulani', 'Gambian', 'Mende',
    'Aari', 'Agaw', 'Amhara', 'Bakola', 'Bedzan', 'Bulala',
    'Elmolo', 'Iraqw', 'Kaba', 'Mada', 'Mursi', 'Ngumba',
    'Ogiek', 'Rendille', 'Sengwer', 'Tikar_South', 'Luo',
    # North African
    'Mozabite', 'Saharawi'
]

# European populations (including groups from Europe, Caucasus, and related)
european_patterns = [
    # Western, Northern, Eastern, Southern European
    'French', 'Sardinian', 'Italian', 'Spanish', 'Basque', 
    'English', 'Scottish', 'Orcadian', 'Russian', 'Adygei', 
    'Czech', 'Croatian', 'Hungarian', 'Lithuanian', 'Estonian',
    'Finnish', 'Norwegian', 'Swedish', 'Icelandic', 'Maltese', 
    'Bulgarian', 'Greek', 'Albanian', 'Tuscan', 'Polish',
    'Bergamo', 'Crete',
    # Caucasian and European-related groups
    'Abkhasian', 'Georgian', 'Chechen', 'North_Ossetian', 'Lezgin',
    'Saami'
]

# Middle Eastern / West Asian (could be considered with European or Asian depending on study)
middle_eastern_patterns = [
    'BedouinB', 'Druze', 'Iraqi_Jew', 'Jordanian', 'Palestinian', 
    'Turkish', 'Yemenite_Jew', 'Iranian', 'Samaritan'
]

# Asian populations (East, South, Southeast, Central, and North Asian)
asian_patterns = [
    # East Asian
    'Han', 'Chinese', 'Dai', 'Japanese', 'Korean', 'Mongolian',
    'Oroqen', 'Hezhen', 'Xibo', 'Miaozu', 'Tujia', 'Yi',
    'Naxi', 'Lahu', 'She', 'Tu',
    # Southeast Asian
    'Cambodian', 'Vietnamese', 'Thai', 'Lao', 'Burmese', 'Dusun', 'Igorot', 'Kinh',
    'Ami', 'Atayal',  # Taiwan aboriginal populations - could be East Asian or Austronesian
    # South Asian
    'Brahmin', 'Kshatriya', 'Vysya', 'Kurumba', 'Irula', 'Bengali', 'Telugu', 'Lambadi', 
    'Chenchu', 'Bonda', 'Santhal', 'Kharia', 'Juang', 'Ho', 'Munda', 'Malayan', 
    'Onge', 'Jarawa', 'Khonda_Dora', 'Mala', 'Punjabi', 'Kapu', 'Madiga', 'Relli', 'Yadava',
    'Kashmiri_Pandit',
    # Central Asian
    'Balochi', 'Brahui', 'Makrani', 'Sindhi', 'Pathan', 'Kalash', 'Burusho', 'Hazara',
    'Tubalar', 'Uygur', 'Kyrgyz', 'Mongola', 'Sherpa', 'Tibetan', 'Tajik', 'Kusunda',
    # North Asian / Siberian populations
    'Yakut', 'Altaian', 'Even', 'Mansi', 'Ulchi', 'Chane'
]

# Northeast Asian / Arctic / Beringian (could be Asian or Indigenous American related)
arctic_patterns = [
    'Chukchi', 'Itelman', 'Eskimo_Chaplin', 'Eskimo_Naukan', 
    'Eskimo_Sireniki', 'Aleut'
]

# Hispanic or Latin American populations (including Indigenous Americas)
latin_american_patterns = [
    # Indigenous Central and South American
    'Maya', 'Pima', 'Colombian', 'Karitiana', 'Surui', 
    'Wichi', 'Piapoco', 'Ticuna', 'Mayan',
    'Mixe', 'Mixtec', 'Nahua', 'Quechua', 'Zapotec',
    # Indigenous North American
    'Chipewyan', 'Cree', 'Tlingit'
]

# Oceanian populations
oceanian_patterns = [
    'Australian', 'Papuan', 'Bougainville', 'Melanesian',
    'Hawaiian', 'Maori'  # Pacific Islander populations
]

# Map populations with a comprehensive approach
if 'Population_ID' in metadata.columns:
    unique_pops = metadata['Population_ID'].unique()
    
    print("\nMapping SGDP populations to ancestry categories:")
    print("-" * 60)
    
    # Initialize counters for category statistics
    category_counts = {
        'African': 0,
        'European': 0,
        'Asian': 0,
        'Hispanic or Latin American': 0,
        'Other/Mixed': 0
    }
    
    # Track populations by category for detailed output
    populations_by_category = {
        'African': [],
        'European': [],
        'Asian': [],
        'Hispanic or Latin American': [],
        'Other/Mixed': []
    }
    
    for pop in unique_pops:
        pop_str = str(pop)
        mapped = False
        
        # Check African patterns
        for pattern in african_patterns:
            if pattern.lower() in pop_str.lower() or pop_str.lower() in pattern.lower():
                sgdp_to_ancestry_category[pop] = 'African'
                print(f"{pop:<25} -> African")
                category_counts['African'] += 1
                populations_by_category['African'].append(pop)
                mapped = True
                break
        
        if not mapped:
            # Check European patterns
            for pattern in european_patterns:
                if pattern.lower() in pop_str.lower() or pop_str.lower() in pattern.lower():
                    sgdp_to_ancestry_category[pop] = 'European'
                    print(f"{pop:<25} -> European")
                    category_counts['European'] += 1
                    populations_by_category['European'].append(pop)
                    mapped = True
                    break
        
        if not mapped:
            # Check Middle Eastern patterns - assign to Asian for consistency with standard practices
            for pattern in middle_eastern_patterns:
                if pattern.lower() in pop_str.lower() or pop_str.lower() in pattern.lower():
                    sgdp_to_ancestry_category[pop] = 'Asian'
                    print(f"{pop:<25} -> Asian (Middle Eastern)")
                    category_counts['Asian'] += 1
                    populations_by_category['Asian'].append(pop)
                    mapped = True
                    break
        
        if not mapped:
            # Check Asian patterns
            for pattern in asian_patterns:
                if pattern.lower() in pop_str.lower() or pop_str.lower() in pattern.lower():
                    sgdp_to_ancestry_category[pop] = 'Asian'
                    print(f"{pop:<25} -> Asian")
                    category_counts['Asian'] += 1
                    populations_by_category['Asian'].append(pop)
                    mapped = True
                    break
        
        if not mapped:
            # Check Arctic patterns - assign to Asian for genetic similarity
            for pattern in arctic_patterns:
                if pattern.lower() in pop_str.lower() or pop_str.lower() in pattern.lower():
                    sgdp_to_ancestry_category[pop] = 'Asian'
                    print(f"{pop:<25} -> Asian (Arctic/Beringian)")
                    category_counts['Asian'] += 1
                    populations_by_category['Asian'].append(pop)
                    mapped = True
                    break
        
        if not mapped:
            # Check Latin American patterns
            for pattern in latin_american_patterns:
                if pattern.lower() in pop_str.lower() or pop_str.lower() in pattern.lower():
                    sgdp_to_ancestry_category[pop] = 'Hispanic or Latin American'
                    print(f"{pop:<25} -> Hispanic or Latin American")
                    category_counts['Hispanic or Latin American'] += 1
                    populations_by_category['Hispanic or Latin American'].append(pop)
                    mapped = True
                    break
        
        if not mapped:
            # Check Oceanian patterns
            for pattern in oceanian_patterns:
                if pattern.lower() in pop_str.lower() or pop_str.lower() in pattern.lower():
                    sgdp_to_ancestry_category[pop] = 'Other/Mixed'
                    print(f"{pop:<25} -> Other/Mixed (Oceanian)")
                    category_counts['Other/Mixed'] += 1
                    populations_by_category['Other/Mixed'].append(pop)
                    mapped = True
                    break
        
        if not mapped:
            # Default to Other/Mixed for truly unmapped populations
            sgdp_to_ancestry_category[pop] = 'Other/Mixed'
            print(f"{pop:<25} -> Other/Mixed (Unknown)")
            category_counts['Other/Mixed'] += 1
            populations_by_category['Other/Mixed'].append(pop)

    print(f"\nTotal populations mapped: {len(sgdp_to_ancestry_category)}")
    print(f"Mapping summary:")
    for category in ['African', 'European', 'Asian', 'Hispanic or Latin American', 'Other/Mixed']:
        count = category_counts[category]
        print(f"  {category}: {count} populations")
    
    # Detailed output of all populations by category
    print("\nDetailed Population Assignments by Ancestry Category:")
    print("="*80)
    for category, pops in populations_by_category.items():
        if pops:
            print(f"\n{category} ({len(pops)} populations):")
            # Print in columns for better readability
            col_width = 25
            cols = 3  # Number of columns for display
            rows = (len(pops) + cols - 1) // cols  # Ceiling division
            
            for i in range(rows):
                row_str = ""
                for j in range(cols):
                    idx = i + j * rows
                    if idx < len(pops):
                        row_str += f"{pops[idx]:<{col_width}}"
                print(row_str)

# --- Load VCF genotype data ---
print("\nLoading SGDP VCF data...")
callset = allel.read_vcf(
    vcf_path,
    fields=['samples', 'calldata/GT', 'variants/POS'],
    log=sys.stdout
)
genotypes = allel.GenotypeArray(callset['calldata/GT'])
positions = callset['variants/POS']
samples = callset['samples']

# --- Extract population from VCF sample names ---
print("Extracting population from VCF sample names...")

def extract_population_from_sample_name(sample_name):
    """
    Extract population name from VCF sample names like 'S_Mozabite-1' or 'B_French-3'
    Returns the population name (e.g., 'Mozabite', 'French')
    """
    # Pattern: [S/B]_PopulationName-Number
    match = re.match(r'[SB]_([^-]+)-\d+', sample_name)
    if match:
        return match.group(1)
    else:
        # Fallback for other patterns
        return sample_name

# Create mapping from VCF samples to populations
vcf_sample_to_population = {}
for sample in samples:
    population = extract_population_from_sample_name(sample)
    vcf_sample_to_population[sample] = population

print(f"Sample extraction examples:")
for i, sample in enumerate(samples[:10]):
    pop = vcf_sample_to_population[sample]
    print(f"  {sample} -> {pop}")

# --- Map VCF samples to ancestry categories ---
print("Mapping VCF samples to ancestry categories...")

sample_to_ancestry_category = []
category_sample_counts = {}
unrecognized_populations = set()

for sample in samples:
    population = vcf_sample_to_population[sample]
    
    # Map population to ancestry category
    if population in sgdp_to_ancestry_category:
        category = sgdp_to_ancestry_category[population]
    else:
        category = 'Other/Mixed'
        unrecognized_populations.add(population)
    
    sample_to_ancestry_category.append(category)
    category_sample_counts[category] = category_sample_counts.get(category, 0) + 1

sample_to_ancestry_category = np.array(sample_to_ancestry_category)

# Print mapping summary
print(f"\nSample mapping summary:")
print(f"Total samples in VCF: {len(samples)}")
print(f"Sample counts by ancestry category:")
for category in sorted(category_sample_counts.keys()):
    print(f"  {category}: {category_sample_counts[category]} samples")

if unrecognized_populations:
    print(f"\nUnrecognized populations (assigned to 'Other/Mixed'):")
    for pop in sorted(unrecognized_populations):
        print(f"  {pop}")

# Get unique categories that have genetic data
unique_categories_with_data = np.unique(sample_to_ancestry_category)
print(f"\nAncestry categories found in genetic data: {unique_categories_with_data}")

# Verify we have samples for major categories
for category in ['African', 'European', 'Asian', 'Hispanic or Latin American']:
    count = np.sum(sample_to_ancestry_category == category)
    print(f"{category}: {count} samples")

# --- Store results ---
pi_mean = {}
pi_ci = {}

# --- Compute nucleotide diversity and confidence intervals ---
print("Computing nucleotide diversity by ancestry category...")

for category in unique_categories_with_data:
    idx = np.where(sample_to_ancestry_category == category)[0]
    if len(idx) == 0:
        continue
        
    print(f"Processing {category}: {len(idx)} samples")
    sub_genotypes = genotypes.take(idx, axis=1)
    
    # Calculate allele counts properly
    ac = sub_genotypes.count_alleles()
    print(f"  Allele count shape: {ac.shape}")
    print(f"  Total variants: {ac.shape[0]}")
    
    # Check for variants with more than 2 alleles (keep only biallelic)
    is_biallelic = ac.max_allele() == 1
    print(f"  Biallelic variants: {np.sum(is_biallelic)}")
    
    if np.sum(is_biallelic) == 0:
        print(f"  Warning: No biallelic variants for {category}")
        continue
    
    # Filter to biallelic sites
    ac_biallelic = ac[is_biallelic]
    pos_biallelic = positions[is_biallelic]
    
    # Check for segregating sites (not monomorphic)
    # For biallelic sites, segregating means both alleles are present
    is_segregating = (ac_biallelic[:, 0] > 0) & (ac_biallelic[:, 1] > 0)
    print(f"  Segregating biallelic variants: {np.sum(is_segregating)}")
    
    if np.sum(is_segregating) == 0:
        print(f"  Warning: No segregating variants for {category}")
        continue
    
    # Use segregating biallelic sites
    ac_final = ac_biallelic[is_segregating]
    pos_final = pos_biallelic[is_segregating]
    
    # Bootstrap settings
    block_size = 5000
    n_boot = 1000

    n_variants = len(pos_final)
    n_blocks = max(1, n_variants // block_size)
    
    print(f"  Using {n_variants} variants in {n_blocks} blocks")
    
    # Calculate blocks
    block_pi = []
    for i in range(n_blocks):
        start_idx = i * block_size
        end_idx = min((i + 1) * block_size, n_variants)
        
        if end_idx > start_idx:
            pos_block = pos_final[start_idx:end_idx]
            ac_block = ac_final[start_idx:end_idx]
            
            if len(pos_block) > 0:
                try:
                    diversity = allel.sequence_diversity(pos_block, ac_block)
                    if not np.isnan(diversity) and diversity > 0:
                        block_pi.append(diversity)
                except Exception as e:
                    print(f"    Error calculating diversity for block {i}: {e}")
                    continue
    
    print(f"  Calculated {len(block_pi)} blocks with valid diversity")
    
    if len(block_pi) == 0:
        print(f"  Warning: No valid diversity blocks for {category}")
        continue

    # Calculate mean diversity
    mean_pi = np.mean(block_pi)
    
    # Bootstrap resampling for confidence intervals
    if len(block_pi) > 1:
        rng = np.random.default_rng(seed=42)
        bootstrap_reps = []
        for _ in range(n_boot):
            sampled = rng.choice(block_pi, size=len(block_pi), replace=True)
            bootstrap_reps.append(np.mean(sampled))
        
        ci_low = np.percentile(bootstrap_reps, 2.5)
        ci_high = np.percentile(bootstrap_reps, 97.5)
        pi_ci[category] = (ci_low, ci_high)
    else:
        # If only one block, use the mean as both bounds
        pi_ci[category] = (mean_pi, mean_pi)
        ci_low, ci_high = mean_pi, mean_pi

    # Store results
    pi_mean[category] = mean_pi
    
    print(f"  π = {mean_pi*10000:.3f} ×10⁻⁴ (95% CI: {ci_low*10000:.3f}-{ci_high*10000:.3f})")

# --- Normalize to get % share of global π ---
total_pi = sum(pi_mean.values())
if total_pi > 0:
    pi_share = {category: pi / total_pi * 100 for category, pi in pi_mean.items()}
else:
    pi_share = {}

print(f"\nNormalized diversity contributions:")
for category, share in pi_share.items():
    print(f"  {category}: {share:.1f}%")

# --- Prepare data for plotting ---
# Define the order of ancestry categories for consistent visualization
ancestry_categories = ['European', 'Asian', 'African', 'African American or Afro-Caribbean', 
                   'Hispanic or Latin American', 'Other/Mixed']

population_pct_values = []
diversity_values = []
category_labels = []

for category in ancestry_categories:
    pop_pct = modern_population_representation.get(category, 0)
    diversity_val = pi_share.get(category, 0)
    
    population_pct_values.append(pop_pct)
    diversity_values.append(diversity_val)
    category_labels.append(category)

# --- Create Figure 1: Population Representation vs Genetic Diversity Contribution ---
fig, ax = plt.subplots(figsize=(16, 10))

x = np.arange(len(category_labels))
width = 0.35

# Create bars
population_bars = ax.bar(x - width/2, population_pct_values, width, label='% of world population', 
                     color='lightblue', alpha=0.8, edgecolor='navy')
diversity_bars = ax.bar(x + width/2, diversity_values, width, label='% of total genetic diversity present in SGDP', 
                        color='orange', alpha=0.8, edgecolor='darkred')

# Customize the plot
ax.set_xlabel('Ancestry Category', fontsize=14, fontweight='bold')
ax.set_ylabel('Percentage of Worldwide Total (%)', fontsize=14, fontweight='bold')
ax.set_title('Disproportionate Contributions to Human Genetic Diversity\nCompared to Modern Global Population Distribution', 
             fontsize=16, fontweight='bold', pad=20)


# Customize x-axis labels with rotation for better readability
ax.set_xticks(x)
ax.set_xticklabels(category_labels, rotation=45, ha='right', fontsize=12)

# Customize legend
ax.legend(fontsize=12, loc='upper right')

# Add grid for better readability
ax.grid(axis='y', alpha=0.3, linestyle='--')
ax.set_axisbelow(True)

# Add value labels on bars
for i, (pop_bar, diversity_bar) in enumerate(zip(population_bars, diversity_bars)):
    # Population percentage labels
    height = pop_bar.get_height()
    if height > 0:
        ax.text(pop_bar.get_x() + pop_bar.get_width()/2., height + 0.5,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Diversity contribution labels
    height = diversity_bar.get_height()
    if height > 0:
        ax.text(diversity_bar.get_x() + diversity_bar.get_width()/2., height + 0.5,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')

# Set y-axis limits
max_val = max(max(population_pct_values), max(diversity_values))
ax.set_ylim(0, max_val * 1.15)

plt.tight_layout()
plt.savefig("Figure1_Disproportionate_Contributions_to_Human_Genetic_Diversity.png", dpi=300, bbox_inches='tight')
plt.show()

# --- Print detailed results ---
print("\n" + "="*100)
print("DISPROPORTIONATE CONTRIBUTIONS TO HUMAN GENETIC DIVERSITY")
print("="*100)
print(f"{'Ancestry Category':<35} {'Modern Global %':<20} {'Genetic Diversity %':<20} {'Contribution Ratio':<18}")
print(f"{'(GWAS classification)':<35} {'(2022 estimates)':<20} {'(Chr 22 π share)':<20} {'(diversity/population)':<18}")
print("-"*100)

for i, category in enumerate(category_labels):
    pop_pct = population_pct_values[i]
    div_pct = diversity_values[i]
    ratio = div_pct / pop_pct if pop_pct > 0 else 0
    print(f"{category:<35} {pop_pct:<20.1f} {div_pct:<20.1f} {ratio:<18.2f}")

# Add explanatory notes with proper framing
print(f"\nEXPLANATION:")
print(f"• Modern Global %: Population distribution based on 2022 UN population estimates")
print(f"• Genetic Diversity %: Proportional contribution to total nucleotide diversity observed in SGDP")
print(f"• Contribution Ratio: Relative genetic diversity contribution compared to global population percentage")
print(f"  - Ratio > 1.0: Population contributes more to genetic diversity than its current global representation")
print(f"  - Ratio < 1.0: Population contributes less to genetic diversity than its current global representation")
print(f"\nMETHODOLOGICAL NOTE: The Simons Genome Diversity Project (SGDP) intentionally samples")
print(f"indigenous and geographically diverse populations to capture human genetic diversity,")
print(f"not to represent current global population distributions. This comparison illustrates how")
print(f"genetic diversity is distributed differently from modern human population density, reflecting")
print(f"evolutionary history rather than recent demographic changes.")
print(f"\nAfrica remains the cradle of human genetic diversity, with indigenous African populations")
print(f"contributing disproportionately to global genetic variation despite representing a")
print(f"smaller percentage of the current world population.")

# --- Print SGDP diversity details ---
print("\n" + "="*90)
print("DETAILED NUCLEOTIDE DIVERSITY ESTIMATES BY ANCESTRY CATEGORY (Chromosome 22)")
print("="*90)
print(f"{'Ancestry Category':<35} {'π (×10⁻⁴)':<15} {'95% CI':<25} {'Share of Total':<15}")
print(f"{'':<35} {'(nucleotide div.)':<15} {'(bootstrap)':<25} {'Diversity %':<15}")
print("-"*90)

for category in ancestry_categories:
    if category in pi_mean:
        pi_val = pi_mean[category] * 10000  # Convert to ×10⁻⁴ scale
        ci_low, ci_high = pi_ci[category]
        ci_low_scaled = ci_low * 10000
        ci_high_scaled = ci_high * 10000
        global_pct = pi_share[category]
        
        print(f"{category:<35} {pi_val:<15.3f} ({ci_low_scaled:.3f}-{ci_high_scaled:.3f}){'':<7} {global_pct:<15.1f}")

# --- Save data ---
results_df = pd.DataFrame({
    'Ancestry_Category': category_labels,
    'Modern_Global_Population_Percent': population_pct_values,
    'Genetic_Diversity_Percent': diversity_values,
    'Contribution_Ratio': [div/pop if pop > 0 else 0 for div, pop in zip(diversity_values, population_pct_values)]
})

results_df.to_csv('Disproportionate_Contributions_to_Human_Genetic_Diversity.csv', index=False)
print(f"\nData saved to 'Disproportionate_Contributions_to_Human_Genetic_Diversity.csv'")

# --- Print population mapping summary ---
print("\n" + "="*60)
print("SGDP SAMPLE DISTRIBUTION BY ANCESTRY CATEGORY")
print("="*60)
category_counts = {}
for sample in samples:
    population = vcf_sample_to_population[sample]
    if population in sgdp_to_ancestry_category:
        category = sgdp_to_ancestry_category[population]
    else:
        category = 'Other/Mixed'
    category_counts[category] = category_counts.get(category, 0) + 1

for category, count in sorted(category_counts.items()):
    print(f"{category}: {count} samples")