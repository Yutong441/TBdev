#' Draw the arrow sign for `rich_forest`
#'
#' @description Please ignore the following arguments: `yval`, `xval`,
#' `group1`, `group2`
#' @param band_ratio how much thicker the arrow should be compared to the
#' `arrow_thickness` field in `AP`.
#' @param length_ratio the ratio of the arrow versus the maximum value along
#' the x axis
#' @param nudge_x by how much the labels for the arrow to shift away from the
#' arrow head.
#' @param nudge_y by how much the arrow and its labels to shift vertically up
#' from the maximum value along the y axis
#' @importFrom ggplot2 aes_string
#' @importFrom magrittr %>%
sagittae <- function (plot_df, aes_param, band_ratio=5, length_ratio=1,
                      nudge_x=0.1, nudge_y=1,
                      yval='yval', xval='emean', 
                      group1='cluster',
                      group2='compare_group'){

        plot_df %>% dplyr::group_by (!!as.symbol (group1)) %>%
                # the arrow ends along the x axis
                dplyr::mutate (finis = pmin (abs (max (!!as.symbol (xval))), 
                                             abs (min (!!as.symbol (xval))) ) ) %>%
                # reset and calculate the y position
                dplyr::ungroup () %>%
                dplyr::group_by (!!as.symbol (group1), !!as.symbol(group2)) %>% 
                # the true length rescaled by `length_ratio`
                dplyr::summarise (max_x = abs(finis)*length_ratio, 
                                  min_x = -abs(finis)*length_ratio,
                                  # adjust the distance of the arrow from the
                                  # main bars
                                  max_y = max (!!as.symbol (yval) ) + nudge_y) -> x_max_df

        # draw the arrow towards the right
        layer1 <- ggplot2::geom_segment(aes_string(x=0, xend='max_x', 
                                                   y='max_y', yend='max_y', 
                                             color=group1), data=x_max_df,
                    arrow = get_arrow (aes_param), #from 'Ax2D.R'
                    size = band_ratio*aes_param$arrow_thickness, 
                    linejoin=aes_param$arrow_linejoin, show.legend=F)

        # draw the arrow towards the left
        layer2 <- ggplot2::geom_segment( aes_string (x=0, xend='min_x', y='max_y', 
                                yend='max_y', color=group2), data=x_max_df,
                    arrow = get_arrow (aes_param), 
                    size = band_ratio*aes_param$arrow_thickness, 
                    linejoin=aes_param$arrow_linejoin, show.legend=F)

        # label the right arrow
        layer3 <- ggplot2::geom_text (aes_string (x='max_x', y='max_y', label=group1, 
                                           color=group1), data=x_max_df,
                vjust='bottom', hjust='left', family=aes_param$font_fam,
                size=aes_param$point_fontsize, nudge_x=nudge_x, show.legend=F)

        # label the left arrow
        layer4 <- ggplot2::geom_text (aes_string (x='min_x', y='max_y', label=group2, 
                                           color=group2), data=x_max_df,
                vjust='bottom', hjust='right', family=aes_param$font_fam,
                size=aes_param$point_fontsize, nudge_x=-nudge_x, show.legend=F)
        return (list(layer1, layer2, layer3, layer4))
}

#' Plot forest of GSEA terms 
#' 
#' @description This function is the backend for `rich_forest`. 
#' @param extend_x_pos how much percentage to extend along the right side of
#' the x axis to prevent the gene label text from being cut off.
#' @param extend_x_neg how much percentage to extend along the left side of
#'
#' the x axis to prevent the gene label text from being cut off.
#' @param ladder_step how much should the two bars by separated from each
#' other. Empirically, I did not find this argument useful.
#' @param band_width the thickness of each bar
#'
#' @param label_shift_ratio by what percentage along the x axis to shift the
#' term labels
#' @param shrink_ratio by how much the fontsize of the labels should decrease
#'
#' @param nudge_y by how much the arrow and its labels to shift vertically up
#' from the maximum value along the y axis
#' @param band_ratio how much thicker the arrow should be compared to the
#' `arrow_thickness` field in `AP`.
#' @param ... arguments for `sagittae`
#' @importFrom ggplot2 aes
#' @importFrom magrittr %>%
gg_silva <- function (plot_df, AP, 
                      # plot
                      extend_x_neg=1., extend_x_pos=1.,
                      # appearance of the bars
                      ladder_step=1, band_width=0.3, 
                      # modify the GSEA/gene labels
                      label_shift_ratio=0.05, shrink_ratio=1., 
                      # arrows
                      nudge_y=1,...){
        plot_df %>% dplyr::group_by (cluster, enriched, compare_group) %>% 
                dplyr::mutate (yval=order (abs (emean) ) ) %>%
                dplyr::mutate (xmax=sign(emean)*max(abs (emean))) %>%
                data.frame () -> plot_df
        
        all_comp <- levels (factor(plot_df$compare_group))
        # the step along the ladder
        plot_df %>% dplyr::count (cluster, enriched, compare_group) %>% 
                dplyr::pull (n) %>% max () -> off_set

        # offset the comparion groups
        for (i in 1:length (all_comp)){
                sel_index <- plot_df$compare_group==all_comp[i]
                plot_df$yval [sel_index] <- plot_df$yval [sel_index] + 
                        (i-1)*(off_set+nudge_y)
        }
        plot_df$yval <- ladder_step*plot_df$yval

        # create a column that contains all the cell types in both group1 and
        # group2 to facilitate coloring
        plot_df$enriched_cell <- as.character (plot_df$compare_group)
        enriched_index <- plot_df$enriched >0
        plot_df$enriched_cell [enriched_index] <- as.character (
                                plot_df$cluster [enriched_index])

        # control the boundary of the plot
        max_x <- (1+ extend_x_pos)*max (plot_df$emean)
        min_x <- (1+ extend_x_neg)*min (plot_df$emean)

        # separate control over the postively and negatively enriched labels
        plot_df$neg_label <- plot_df$glabel
        plot_df$glabel [!enriched_index] <- NA
        plot_df$neg_label [enriched_index] <- NA

        ggplot2::ggplot (plot_df, aes (x=emean, y= yval)) +
                ggplot2::geom_bar (aes (fill=enriched_cell), stat='identity',
                                   show.legend=F, width=band_width, orientation='y')+ 

                # add text for the tree branches on the right
                ggplot2::geom_text (aes(label= glabel, 
                                        x=emean+xmax*label_shift_ratio), 
                                    family=AP$font_fam, hjust='left', 
                                    size=AP$point_fontsize*shrink_ratio)+

                # add text for the tree branches on the left
                ggplot2::geom_text (aes(label=neg_label, 
                                        x=emean+xmax*label_shift_ratio), 
                                    family=AP$font_fam, hjust='right', 
                                    size=AP$point_fontsize*shrink_ratio)+

                # add arrow sign
                sagittae (plot_df, AP, nudge_y=nudge_y,...)+
                ggplot2::facet_wrap (~cluster, scales='free_y')+

                # central dashed line
                ggplot2::geom_vline (xintercept=0, linetype='dashed', 
                                     size=AP$pointsize/3, color='black') +

                # apply style 
                theme_TB ('dotplot', feature_vec=plot_df$enriched_cell, rotation=0,
                          color_fill=T, aes_param=AP)+
                ggplot2::theme (axis.text.y=ggplot2::element_blank (), 
                                strip.text=ggplot2::element_blank ()) +
                add_custom_color (feature_vec=c(as.character (plot_df$cluster),
                                                as.character (plot_df$compare_group)),
                                  color_fill=F, aes_param=AP)+
                custom_tick (x_y='x')+

                # coordinates
                ggplot2::ylab('')+ ggplot2::labs(fill='p value') +
                ggplot2::xlab ('enrichment score') +
                ggplot2::xlim (c(min_x, max_x))+
                # prevent labels from being cut off in facets
                ggplot2::coord_cartesian (clip='off') 
}

#' Barplot for GSEA results
#' 
#' @description This function creates a compact way of visualising GSEA
#' results by creating multiple tree-like plots, or 'forests' (silva in Latin).
#' Each tree represents a pairwise comparison between one particular cell type
#' with all other cell types. This pattern is repeated for other cell types as
#' well. Each cell pair comparison is marked by an arrow (sagitta in Latin) on
#' top 
#' @param plot_data a dataframe generated from `run_GSEA_all_types`. NB: if the
#' input data frame is from other functions, make sure you apply some p value
#' cut off. It must have the following columns:
#' @param category_col the GSEA terms
#' @param enrich_col the enrichment value
#' @param pval_col the (adjusted) p value
#' @param term_col the ID for each gene in the format of A.B.entrez where A and
#' B can be of anything and entrez means entrez ID
#' @param group1 comparison group1. Each group in `cluster` will be in a
#' @param separate tree. You may change it by specifying `group`=...`
#' @param group2 comparison group2. Each group here will
#' be in the same tree.
#'
#' @param markers a dataframe generated from `find_DE_genes` that contain the
#' log fold change of gene expression. This information is used to order the
#' appearance of genes in the plot labels. It must contain the columns:
#' `group` (group1), `logFC` (used to rate which genes to show), `feature`
#' (genes), and `compare_group` (group2). You can change `compare_group` via:
#' @param compare_group_col_marker which column in `markers` contain the group2
#' information.
#' @param organism_db a gene database
#'
#' @param show_num how many terms to show per pairwise comparison
#' @param show_gene_labels how many genes to show after each term
#'
#' @param simplification whether to simplify the GO/KEGG/Reactome terms
#' @param sim_dict a data frame with 2 columns: 'ori' for the original terms,
#' 'sub' for the strings that will replace the original terms. The default is a
#' a limited data frame built into this package.
#' @param append_default_dict append default dictionary to simplify the terms
#' not exist.
#' @param AP aesthetic parameters for plotting
#' @param ... other parameters to pass onto `gg_silva` and `sagittae`. 
#' @export
rich_forest <- function (plot_data, markers, organism_db, 
                         # genes and labels
                         show_num=4,
                         show_gene_labels=3, 

                         # columns in `plot_data`
                         category_col='category',
                         pval_col='p.adjust',
                         enrich_col='value',
                         term_col='termID',
                         group1='cluster',
                         group2='compare_group',

                         # simplify and clean GSEA terms
                         simplification=T,
                         sim_dict=NULL, append_default_dict=T,

                         compare_group_col_marker='compare_group',
                         AP=NULL, ...){
        AP <- return_aes_param (AP)
        sim_dict <- append_default_dictionary (sim_dict, append_default_dict)
        relevant_terms <- remove_terms (plot_data [, category_col], AP)
        plot_data %>% dplyr::filter (!!as.symbol (category_col) %in% 
                                     relevant_terms) -> plot_data
        if (simplification) {plot_data <- simplify_gsea (plot_data, sim_dict, 
                term_col=category_col)} # from 'clean_terms.R'

        # make sure the dataframe has a column called 'compare_group'
        markers$compare_group <- markers [, compare_group_col_marker]
        plot_data$compare_group <- plot_data[, group2]

        # For comparison of one cell type with the rest of the cell types,
        # clean up the GSEA terms and append the gene names
        plot_data %>% dplyr::pull (compare_group) %>% unique () -> all_comp
        pdata_list <- list ()
        for (i in 1:length(all_comp)){
                sel_plot <- plot_data [plot_data$compare_group==all_comp[i],]
                sel_mark <- markers [markers$compare_group==all_comp[i],]
                # pick the most significant GSEA terms. See `enrichment.R`
                plot_df <- summarise_gsea (sel_plot, show_num, category_col=category_col, 
                                           pval_col=pval_col, enrich_col=enrich_col, 
                                           cluster_col=group1)
                # append the genes behind the GSEA terms and take care of the
                # line breaks. See `line_break.R`
                plot_df$glabel <- term_gene_labels (sel_plot, plot_df, organism_db,
                                show_gene_labels, markers, category_col=category_col, 
                                cluster_col=group1, term_col=term_col)
                pdata_list[[i]] <- plot_df
        }
        plot_df <- do.call (rbind, pdata_list)
        # ensure the groups appear in the desired order
        plot_df [, group1]<- partial_relevel (plot_df [, group1], AP$cell_order)
        group2_lev <- partial_relevel (plot_df [, group2], AP$cell_order)
        plot_df [, group2] <- factor (group2_lev, levels=rev (levels (group2_lev)) )
        return (gg_silva (plot_df, AP, group1=group1,...))
}
