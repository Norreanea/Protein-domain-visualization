# Load necessary libraries
library(Biostrings)
library(msa)
library(seqinr)
library(dplyr)
library(drawProteins)
library(ggplot2)
library(ggpubr)
library(ggmsa)
library(cowplot)
library(httr)
library(UniProt.ws)
library(ggrepel)
library(gridExtra)
#install.packages("ggthemes")
library(ggthemes)
library(scales)
library(plyr)
library(tidyr)
library(purrr)
# Function to read sequences and subset them based on unique names
read_and_subset_sequences = function(fname) {
  file = readAAStringSet(fname)
  unique_RNase = unique(sub("\\s+\\w+\\s+\\w+$", "", names(file)))
  sequences = list()
  for (i in seq_along(unique_RNase)) {
    seq_name = unique_RNase[i]
    sequences[[seq_name]] = file[names(file)[startsWith(names(file), seq_name)]]
  }
  return(sequences)
}

# Function to adjust sequence names
adjust_sequence_names = function(seq, start_pos) {
  names(seq) = substr(names(seq), start_pos, nchar(names(seq)))
  names(seq) = gsub(" ", "_", names(seq))
  return(seq)
}

# Function to perform MSA and calculate distance
perform_msa_and_distance = function(seq) {
  aln1 = msa(seq, "ClustalOmega")
  aln2 = msaConvert(aln1, type = "seqinr::alignment")
  dist = dist.alignment(aln2, "identity")
  return(list(alignment = aln1, distance = dist))
}

# Function to map amino acid positions
mapAminoAcidPositions = function(original_seq, aligned_seq) {
  original_aas = unlist(strsplit(original_seq, ""))
  aligned_aas = unlist(strsplit(aligned_seq, ""))
  original_positions = numeric(0)
  aligned_positions = numeric(0)
  original_seq_amino = character(0)
  aligned_seq_amino = character(0)
  original_pos_counter = 1
  
  for (i in seq_along(aligned_aas)) {
    if (aligned_aas[i] != "-") {
      if (original_pos_counter <= length(original_aas)) {
        original_positions = c(original_positions, original_pos_counter)
        original_seq_amino = c(original_seq_amino, original_aas[original_pos_counter])
        aligned_seq_amino = c(aligned_seq_amino, aligned_aas[i])
        aligned_positions = c(aligned_positions, i)
        original_pos_counter = original_pos_counter + 1
      }
    }
  }
  
  result = data.frame(
    Original_seq_amino = original_seq_amino,
    Original_position = original_positions,
    Aligned_seq_amino = aligned_seq_amino,
    Aligned_position = aligned_positions
  )
  return(result)
}

# Function to fetch domain data from PFAM API
fetch_pfam_domains = function(uniprot_ids) {
  base_url = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/UniProt/"
  ids = strsplit(uniprot_ids, " ")[[1]]
  
  domain_list = list()
  gene_list = list()
  for (id in ids) {
    url = paste0(base_url, id, "/?page_size=100")
    response = httr::GET(url)
    data = httr::content(response, "parsed")
    uniprot_gene = UniProt.ws::mapUniProt("UniProtKB_AC-ID", "UniProtKB", query = id)
    
    split_names = strsplit(uniprot_gene$Gene.Names, " ")[[1]]
    first_gene_name = split_names[1]
    
    gene_list[[id]] = data.frame(
      gene = first_gene_name,
      accession = id,
      stringsAsFactors = FALSE
    )
    
    if (!is.null(data$results)) {
      for (result in data$results) {
        for (protein in result$proteins) {
          for (location in protein$entry_protein_locations) {
            for (fragment in location$fragments) {
              domain_list[[length(domain_list) + 1]] = list(
                type = "DOMAIN",
                description = result$metadata$name,
                begin = as.numeric(fragment$start),
                end = as.numeric(fragment$end),
                accession = id,
                taxid = as.integer(protein$organism)
              )
            }
          }
        }
      }
    }
  }
  
  gene_df = do.call(rbind, gene_list)
  
  if (length(domain_list) > 0) {
    domain_df = do.call(rbind, lapply(domain_list, function(x) data.frame(x, stringsAsFactors = FALSE)))
  } else {
    domain_df = data.frame(type = character(), description = character(), begin = integer(), end = integer(), accession = character(), taxid = integer(), stringsAsFactors = FALSE)
  }
  
  return(list(domain_df = domain_df, gene_df = gene_df))
}

# Function to add non-overlapping domains from PFAM API into prot_data
add_non_overlapping_domains = function(prot_data, new_domain_data) {
  non_overlapping_domains = new_domain_data[!sapply(seq(nrow(new_domain_data)), function(i) {
    any(prot_data$begin <= new_domain_data$end[i] & prot_data$end >= new_domain_data$begin[i] & prot_data$type == "DOMAIN")
  }), ]
  
  combined_data = dplyr::bind_rows(prot_data, non_overlapping_domains)
  combined_data = combined_data[order(combined_data$begin), ]
  return(combined_data)
}

# Function to visualize protein domains for given UniProt IDs
collect_protein_domains = function(rnase_name, uniprot_ids, order, title, msa_results) {
  ids = unlist(strsplit(uniprot_ids, " "))
  prot_data = drawProteins::get_features(uniprot_ids)
  prot_data = drawProteins::feature_to_dataframe(prot_data)
  new_domain_data = fetch_pfam_domains(uniprot_ids)
  prot_data = add_non_overlapping_domains(prot_data, new_domain_data$domain_df)
  
  taxid_to_species = list(
    "9606" = "Homo_sapiens",
    "10090" = "Mus_musculus",
    "559292" = "Saccharomyces_cerevisiae",
    "6239" = "Caenorhabditis_elegans",
    "7227" = "Drosophila_melanogaster"
  )
  prot_data$scientific_name = sapply(as.character(prot_data$taxid), function(x) taxid_to_species[[x]])
  
  scientific_to_short = list(
    "Homo_sapiens" = "HS",
    "Mus_musculus" = "MM",
    "Saccharomyces_cerevisiae" = "SC",
    "Caenorhabditis_elegans" = "CE",
    "Drosophila_melanogaster" = "DM"
  )
  prot_data$short_tax = sapply(prot_data$scientific_name, function(x) scientific_to_short[[x]])
  #prot_data$entryName = paste(prot_data$short_tax, prot_data$gene, sep = "_")
  if (rnase_name != "AGO2") {
    original_seq = msa_results[[rnase_name]]$alignment@unmasked
    aligned_seq = msa_results[[rnase_name]]$alignment
    for (i in seq_along(prot_data$accession)) {
      species_name = taxid_to_species[[as.character(prot_data$taxid[i])]]
      if (!is.null(species_name)) {
        position_mapping = mapAminoAcidPositions(as.character(original_seq[[species_name]]), as.character(aligned_seq@unmasked[[species_name]]))
        prot_data$begin[i] = position_mapping$Aligned_position[match(prot_data$begin[i], position_mapping$Original_position)]
        prot_data$end[i] = position_mapping$Aligned_position[match(prot_data$end[i], position_mapping$Original_position)]
      }
    }
  } else {
    original_seq = msa_results[["AGO2_full_raw"]]
    aligned_seq = msa_results[["AGO2Aln1_full"]]
    entryName_to_species = list(
      "G5EES3" = "CE_alg-1",
      "O16720" = "CE_alg-2",
      "Q9UKV8" = "HS_AGO2",
      "Q9UL18" = "HS_AGO1"
    )
    for (i in seq_along(prot_data$accession)) {
      species_name = entryName_to_species[[as.character(prot_data$accession[i])]]
      if (!is.null(species_name)) {
        position_mapping = mapAminoAcidPositions(as.character(original_seq[[species_name]]), as.character(aligned_seq@unmasked[[species_name]]))
        prot_data$begin[i] = position_mapping$Aligned_position[match(prot_data$begin[i], position_mapping$Original_position)]
        prot_data$end[i] = position_mapping$Aligned_position[match(prot_data$end[i], position_mapping$Original_position)]
      }
    }
  }
  
  if (length(order) == length(ids) && all(order <= length(ids))) {
    prot_data$order = order[match(prot_data$accession, ids)]
  } else {
    cat("Order vector length does not match or invalid order indices")
  }
  prot_data = merge(prot_data, new_domain_data$gene_df, by = "accession", all = TRUE)
  prot_data$entryName = paste(prot_data$short_tax, prot_data$gene, sep = "_")
  prot_data$set = title

  
  return(prot_data)
}

# Function to add labels to a plot
add_label = function(plot, label, xpos, ypos) {
  plot + annotate("text", label = label, x = xpos, y = ypos, hjust = 0.5, vjust = 0.5, 
                  color = "black", size = 5, fontface = "bold")
}

# Main 
fname = "D:/Review/all_RNases_new.fa"
sequences = read_and_subset_sequences(fname)

# Adjust sequence names
start_positions = c(17, 27, 27, 27, 38, 26, 23, 17, 33, 31, 40)
for (i in seq_along(sequences)) {
  sequences[[i]] = adjust_sequence_names(sequences[[i]], start_positions[i])
}

# Perform MSA and calculate distance
msa_results = list()
for (seq_name in names(sequences)) {
  msa_results[[seq_name]] = perform_msa_and_distance(sequences[[seq_name]])
}
#msa_results$`argonaute RISC catalytic component 2`
# Exception:
AGO2_full="D:/Review/AGO2.fa"
AGO2_full_raw=readAAStringSet(AGO2_full)
AGO2Aln1_full=msa(readAAStringSet(AGO2_full),"ClustalOmega")
msa_results$AGO2$alignment=AGO2Aln1_full
msa_results$AGO2_full_raw=AGO2_full_raw
msa_results$AGO2Aln1_full=AGO2Aln1_full
names(msa_results)

# Collect information
RNASEH2A_dom = collect_protein_domains("ribonuclease H2 subunit A", "O75792 Q9CWY8 P53942", c(3, 2, 1), "Ribonuclease H2 subunit A", msa_results)
RNASEH2B_dom = collect_protein_domains("ribonuclease H2 subunit B", "Q5TBB1 Q80ZV0", c(2, 1), "Ribonuclease H2 subunit B", msa_results)
AGO2_dom = collect_protein_domains("AGO2", "G5EES3 O16720 Q9UL18 Q9UKV8", c(2, 1, 4, 3), "Argonaute-2", msa_results)
DICER_dom = collect_protein_domains("dicer 1 ribonuclease III","Q9UPY3 Q9VCU9", c(2, 1), "Endoribonuclease Dicer", msa_results)
ELAC2_dom = collect_protein_domains("elaC ribonuclease Z 2","Q9BQ52 P36159 Q8MKW7 Q80Y81", c(4,1,2, 3), "Zinc phosphodiesterase ELAC protein 2", msa_results)
DIS3L2_dom = collect_protein_domains("DIS3 like 3-5 exoribonuclease 2","Q8IYB7 Q8CI75", c(2, 1), "DIS3-like exonuclease 2", msa_results)

#AGO2_dom[AGO2_dom$type == "CHAIN",]
domains_merged = rbind(RNASEH2A_dom, RNASEH2B_dom, AGO2_dom, DICER_dom, ELAC2_dom, DIS3L2_dom)


#prot_data = domains_merged
# Reorder
desired_order = c("Ribonuclease H2 subunit A", "Ribonuclease H2 subunit B", "Argonaute-2", 
                   "Endoribonuclease Dicer", "DIS3-like exonuclease 2", "Zinc phosphodiesterase ELAC protein 2")
domains_merged$set = factor(domains_merged$set, levels = desired_order)

# Add selected mutation positions
mutation_pos = data.frame(set = c("Ribonuclease H2 subunit A","Ribonuclease H2 subunit B","Ribonuclease H2 subunit B",
                                  "Argonaute-2", "Endoribonuclease Dicer", "Endoribonuclease Dicer", "DIS3-like exonuclease 2",
                                  "Zinc phosphodiesterase ELAC protein 2", "Zinc phosphodiesterase ELAC protein 2", "Zinc phosphodiesterase ELAC protein 2"),
                          begin = c(93,203,234,318,1205,1206,707,198,678,700),
                          order = c(3,2,2,4,2,2,2,4,4,4),
                          type = rep("MUT",10))

domains_merged_with_mutations = rbind.fill(domains_merged, mutation_pos)
domains_merged_with_mutations[domains_merged_with_mutations$type=="MUT",]

# unique(domains_merged$set)
# domains_merged
# "The following columns were added to unify the starting points for CHAIN labels:
max_begin_chain=aggregate(domains_merged_with_mutations$begin[domains_merged_with_mutations$type=="CHAIN"] ~ 
                            domains_merged_with_mutations$set[domains_merged_with_mutations$type=="CHAIN"], FUN = "max")
colnames(max_begin_chain)=c("set","max_begin_chain")
max_end_chain=aggregate(domains_merged_with_mutations$end[domains_merged_with_mutations$type=="CHAIN"] ~ 
                            domains_merged_with_mutations$set[domains_merged_with_mutations$type=="CHAIN"], FUN = "max")
colnames(max_end_chain)=c("set","max_end_chain")



domains_merged_with_mutations = merge(domains_merged_with_mutations,max_begin_chain, by = "set", all = TRUE)
domains_merged_with_mutations = merge(domains_merged_with_mutations,max_end_chain, by = "set", all = TRUE)
domains_merged_with_mutations$description[domains_merged_with_mutations$description=="Ydr279p protein family (RNase H2 complex component) wHTH domain"] =
  "Ydr279p protein family (RNase H2 complex component) \nwHTH domain"
domains_merged_with_mutations[domains_merged_with_mutations$type =="DOMAIN",]
# Rename PAZ domains for uniqueness
domains_merged_with_mutations = domains_merged_with_mutations %>%
  mutate(description = case_when(
    set == "Argonaute-2" & description == "PAZ" ~ "PAZ ",
    TRUE ~ description
  ))
                               
# Custom colors for domains
my_darker_colors = c(
  'N-terminal domain of argonaute' = "#194d77",
  'PAZ' = "#cc660b",
  'PAZ ' = "#cc660b",
  'Argonaute linker 1 domain' = "#238022",
  'Piwi' = "#ac2020",
  'RNB domain' = "#764f99",
  'Rrp44-like cold shock domain' = "#704b3d",
  'DIS3-like exonuclease 2 C terminal' = "#b55e9b",
  'Dis3-like cold-shock domain 2 (CSD2)' = "#7fcccc",
  'Helicase ATP-binding' = "#97901a",
  'Dicer, partner-binding domain' = "#1296a1",
  'Helicase C-terminal' = "#8ba4ba",
  'Dicer dsRNA-binding fold' = "#cc975e",
  'RNase III 1' = "#79b274",
  'RNase III 2' = "#cc7a78",
  'DRBM' = "#9e8ca7",
  'RNase H type-2' = "#9d7d77",
  'Ydr279p protein triple barrel domain' = "#c596a9",
  'Ydr279p protein family (RNase H2 complex component) \nwHTH domain' = "#cc527a",
  'tRNase Z endonuclease' = "#b0b073",
  'Beta-lactamase superfamily domain' = "#7bb0b7"
)

my_darker_colors_df=as.data.frame(my_darker_colors)
colnames(my_darker_colors_df)="color"
my_darker_colors_df$description = rownames(my_darker_colors_df)
domains_merged_with_mutations = merge(domains_merged_with_mutations,my_darker_colors_df, by = "description", all = TRUE)

desired_order = c("Ribonuclease H2 subunit A", "Ribonuclease H2 subunit B", "Argonaute-2", 
                  "Endoribonuclease Dicer", "DIS3-like exonuclease 2", "Zinc phosphodiesterase ELAC protein 2")
domains_merged_with_mutations$set = factor(domains_merged_with_mutations$set, levels = desired_order)

# Define custom levels for each set
custom_levels_list = list(
  "Ribonuclease H2 subunit A" = c("RNase H type-2"),
  "Ribonuclease H2 subunit B" = c("Ydr279p protein triple barrel domain",
                                  "Ydr279p protein family (RNase H2 complex component) \nwHTH domain"),
  "Argonaute-2" = c("N-terminal domain of argonaute", "Argonaute linker 1 domain",
                    "PAZ ", "Piwi"),
  "Endoribonuclease Dicer" = c("Helicase ATP-binding", "Dicer, partner-binding domain",
                               "Helicase C-terminal", "Dicer dsRNA-binding fold",
                               "PAZ", "RNase III 1", "RNase III 2", "DRBM"),
  "DIS3-like exonuclease 2" = c("Rrp44-like cold shock domain", "Dis3-like cold-shock domain 2 (CSD2)",
                                "RNB domain", "DIS3-like exonuclease 2 C terminal"),
  "Zinc phosphodiesterase ELAC protein 2" = c("tRNase Z endonuclease", "Beta-lactamase superfamily domain")
)

# Define a function to set factor levels for a subset of data
set_custom_levels = function(df, custom_levels_list) {
  set_name <- unique(df$set)
  if (length(set_name) == 1 && set_name %in% names(custom_levels_list)) {
    df$description <- factor(df$description, levels = custom_levels_list[[set_name]])
  }
  return(df)
}

# Apply the custom function to each subset and recombine the data
domains_merged_with_mutations = domains_merged_with_mutations %>%
  group_split(set) %>%
  map_dfr(~ set_custom_levels(.x, custom_levels_list))
                               
final_plot_protein_domains_with_legend = by(data = domains_merged_with_mutations, INDICES = domains_merged_with_mutations$set, FUN = function(m) {
  
  m = droplevels(m)
  m = ggplot2::ggplot() +
    ggplot2::ylim(0.5, 6) +
    ggplot2::xlim(-max(m$end, na.rm = TRUE) * 0.4,
                  max(m$end, na.rm = TRUE) + max(m$end, na.rm = TRUE) * 0.1) +
    ggplot2::labs(x = "Amino acid number") +
    ggplot2::labs(y = "") + ggplot2::labs(fill="") +
    guides(fill=guide_legend(nrow =4)) +
    theme_classic(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          axis.ticks = element_blank(), axis.text.y = element_blank(),
          panel.border = element_blank(),
          legend.position = "bottom") +  # Adjust plot margins
    scale_x_continuous(labels = function(x) ifelse(x < 0, "", x)) +  # Hide negative labels
    #facet_wrap(~ set, scales = "free", strip.position = "left", ncol = 2, shrink = TRUE) +
    ggplot2::geom_rect(data = m[m$type == "CHAIN", ],
                       mapping = ggplot2::aes(xmin = begin,
                                              xmax = end,
                                              ymin = order - 0.2,
                                              ymax = order + 0.2),
                       colour = "black",
                       fill = "grey",
                       size = 0.5,
                       alpha = 1.0) +
    ggplot2::geom_text(data = m[m$type == "CHAIN", ],
                       mapping = ggplot2::aes(x = max_begin_chain - (max_end_chain * 0.3),  # Adjust x position to ensure labels are visible
                                              y = order,
                                              label = entryName),
                       hjust = 0.1,
                       vjust = 0.5,
                       size = 3) +
    ggplot2::geom_rect(data = m[m$type == "DOMAIN", ],
                       mapping = ggplot2::aes(xmin = begin,
                                              xmax = end,
                                              ymin = order - 0.25,
                                              ymax = order + 0.25,
                                              fill = description),
                       alpha = 1.0,
                       show.legend = TRUE) +
    scale_fill_manual(values = my_darker_colors) +
    geom_point(data = m[m$type == "MUT", ], 
               mapping = aes(x = begin, y = order + 0.4), size = 3, shape = 25, fill = "red")+
    theme(
      strip.background = element_blank(),
      strip.text = element_blank()
    )
  
  #return(p)
})


final_plot = do.call(cowplot::plot_grid, c(final_plot_protein_domains_with_legend, ncol = 2, align = "hv",labels = "AUTO"
                              ))

# Define the necessary data for each RNase
rnases_data = list(
  list(name = "ribonuclease H2 subunit A", organisms = c("Homo_sapiens", "Mus_musculus", "Saccharomyces_cerevisiae"), new_names = c("HS", "MM", "SC"), start = 90, end = 95, highlight = 93),
  list(name = "ribonuclease H2 subunit B", organisms = c("Homo_sapiens", "Mus_musculus"), new_names = c("HS", "MM"), start = c(198, 234), end = c(203, 239), highlight = c(203, 234)),
  list(name = "AGO2", organisms = c("HS_AGO1", "HS_AGO2", "CE_alg-1", "CE_alg-2"), new_names = c("HS_AGO1", "HS_AGO2", "CE_alg-1", "CE_alg-2"), start = 315, end = 320, highlight = 318),
  list(name = "dicer 1 ribonuclease III", organisms = c("Homo_sapiens", "Drosophila_melanogaster"), new_names = c("HS", "DM"), start = 1205, end = 1210, highlight = c(1205, 1206)),
  list(name = "DIS3 like 3-5 exoribonuclease 2", organisms = c("Homo_sapiens", "Mus_musculus"), new_names = c("HS", "MM"), start = 707, end = 767, highlight = c(707:767)),
  list(name = "elaC ribonuclease Z 2", organisms = c("Homo_sapiens", "Mus_musculus", "Drosophila_melanogaster", "Saccharomyces_cerevisiae"), new_names = c("HS", "MM", "DM", "SC"), start = c(195, 675, 700), end = c(200, 680, 705), highlight = c(198, 678, 700))
)

create_ggmsa_plot = function(alignment, start, end, highlight, rnase_name ) {
  plot_list <- list()
  for (i in seq_along(start)) {
    if(rnase_name %in% c("DIS3 like 3-5 exoribonuclease 2","dicer 1 ribonuclease III")){
      highlight_range <- highlight
      plot <- ggmsa(alignment, color = "LETTER", seq_name = TRUE, by_conservation = TRUE,
                    start = start[i], end = end[i], char_width = 0.9,
                    font = "TimesNewRoman", position_highlight = highlight_range) +
        theme(text = element_text(size = 10),
              plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
        scale_x_continuous(breaks = highlight_range[i])
      plot_list[[i]] <- plot
    } else {
      highlight_range <- highlight
      plot <- ggmsa(alignment, color = "LETTER", seq_name = TRUE, by_conservation = TRUE,
                    start = start[i], end = end[i], char_width = 0.9,
                    font = "TimesNewRoman", position_highlight = highlight_range) +
        theme(text = element_text(size = 10),
              plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
        scale_x_continuous(breaks = highlight_range)
      plot_list[[i]] <- plot
    }
    
  }
  return(plot_list)
}

# Function to generate alignment plots for a given RNase
generate_alignment_plots = function(rnase_name, msa_results) {
  rnase_data = lapply(rnases_data, function(x) if (x$name == rnase_name) x else NULL)
  rnase_data = rnase_data[!sapply(rnase_data, is.null)][[1]]
  
  if (is.null(rnase_data)) {
    stop(paste("RNase", rnase_name, "not found in the provided data."))
  }
  
  aln_str = as(msa_results[[rnase_name]]$alignment, "XStringSet")
  order_index = na.omit(match(rnase_data$organisms, names(aln_str)))
  selected_aln = aln_str[order_index]
  
  if (length(order_index) > 0) {
    names(selected_aln) = rnase_data$new_names[seq_along(order_index)]  # Ensure names are properly assigned and unique
    plots = create_ggmsa_plot(selected_aln, rnase_data$start, rnase_data$end, rnase_data$highlight, rnase_name)
    return(plots)
  } else {
    warning(paste("No matching sequences found for", rnase_name))
    return(NULL)
  }
}


selected_RNASEH2A_aln = generate_alignment_plots("ribonuclease H2 subunit A", msa_results)
selected_RNASEH2B_aln = generate_alignment_plots("ribonuclease H2 subunit B", msa_results)
selected_AGO2_aln = generate_alignment_plots("AGO2", msa_results)
selected_DICER1_aln = generate_alignment_plots("dicer 1 ribonuclease III", msa_results)
selected_DIS3L2_aln = generate_alignment_plots("DIS3 like 3-5 exoribonuclease 2", msa_results)
selected_ELAC2_aln = generate_alignment_plots("elaC ribonuclease Z 2", msa_results)


selected_align = ggdraw() + draw_plot(selected_RNASEH2A_aln[[1]], x = 0.11, y = .91, width = .13, height = .07) +
 draw_plot(selected_RNASEH2B_aln[[1]], x = 0.7, y = .877, width = .13, height = .07) +
 draw_plot(selected_RNASEH2B_aln[[2]], x = 0.82, y = .877, width = .13, height = .07) +
 draw_plot(selected_AGO2_aln[[1]], x = 0.1, y = .605, width = .17, height = .08) +
 draw_plot(selected_DICER1_aln[[1]], x = 0.76, y = .544, width = .13, height = .07) +
 draw_plot(selected_ELAC2_aln[[1]], x = 0.57, y = .27, width = .14, height = .065) +
 draw_plot(selected_ELAC2_aln[[2]], x = 0.71, y = .27, width = .14, height = .065) +
 draw_plot(selected_ELAC2_aln[[3]], x = 0.85, y = .27, width = .14, height = .065) +
 draw_plot(selected_DIS3L2_aln[[1]], x = 0.03, y = -0.26, width = .5, height = 1) 

# Combine the final plot with alignment plots

(combined_plot = ggdraw() +
  draw_plot(final_plot, 0, 0, 1, 1) +
  draw_plot(selected_align, 0, 0, 1, 1))
