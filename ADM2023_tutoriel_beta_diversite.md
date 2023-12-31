R Notebook
================

``` bash
git config --global user.email "you@example.com"
git config --global user.name "Your Name"
```

1.  On commence par charger toutes les librairies nécessaires. On charge
    également tous les outils en indiquant à R le chemin d’accès.

``` r
library(phyloseq)
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
devtools::load_all(path="/home/rstudio/ADM2023_tutoriel/course-material-main/R")
```

    ## ℹ Loading ANF_metaB

On commence par créer une variable output_beta qui contient le chemin
d’accès vers un fichier. On demande à R de créer une direction vers la
variable output_beta portant ce nom si jamais il n’existe pas.
L’argument recursive = TRUE signifie qu’on autorise le fichier
output_beta à créer des sous-dossiers.

``` r
output_beta <- here::here("outputs", "beta_diversity")
if (!dir.exists(output_beta)) dir.create(output_beta, recursive = TRUE)
```

On copie le dossier asv-table contenu dans le course-material-main dans
le dossier data contenu dans ADM2023_tutoriel.

``` bash
cp -R course-material-main/data/asv_table ./data/
```

On crée une variable physeq à laquelle on attribue la fonction pour lire
le fichier .rds.

``` r
physeq <- readRDS(here::here("data",
                             "asv_table",
                             "phyloseq_object_alpha_beta_div.rds"))
```

2.  

On a deux approches différentes pour le processus de normalisation des
librairies. La première consiste à faire un échantillonage pour réduire
le nombre d’observation. Cependant, cela revient à réduire
l’échantillonage donc cela faire perdre des données. La seconde méthode
emploie une analyse de données avec des ratios logarythmiques.

Première méthode : On regarde combien de séquences on a par échantillon
et on les traduit sous forme de tableau pour regarder leur abondance
relative.

``` r
rowSums(physeq@otu_table@.Data)
```

    ## S11B  S1B  S2B  S2S  S3B  S3S  S4B  S4S  S5B  S5S  S6B  S6S  S7B  S7S  S8B  S8S 
    ##  975  837  893  983  878  889  917 1077 1018 1006 1076  937  878  936  846  958 
    ##  S9B  S9S 
    ##  888  991

On crée un objet readsumsdf à qui on attribue la fonction data.frame,
cela fait un tableau. nreads : nombre de séquences lues. decreasing =
TRUE : les séquences sont classées par ordre décroissant d’abondance.
sort = 1 : les séquences sont numérotées de 1 au nombre de taxa de
physed. type = OTU : met les séquences dans des OTU.

Ensuite, on crée un objet tmp à qui on attribue la fonction data.frame.
On applique le même genre de paramètres pour créer un tableau similaire
à readsumsdf.

On combine les tableaux obtenus en un seul tableau avec la fonction
rbind. Puis on affiche uniquement les séquences les plus abondantes avec
la fonction head().

``` r
readsumsdf <- data.frame(nreads = sort(taxa_sums(physeq), decreasing = TRUE),
                        sorted = 1:ntaxa(physeq),
                        type = "OTUs")

tmp <- data.frame(nreads = sort(sample_sums(physeq), decreasing = TRUE), 
                  sorted = 1:nsamples(physeq),
                  type = "Samples")

readsumsdf <- rbind(readsumsdf, tmp)

head(readsumsdf)
```

    ##      nreads sorted type
    ## ASV1   1558      1 OTUs
    ## ASV2    973      2 OTUs
    ## ASV3    899      3 OTUs
    ## ASV4    833      4 OTUs
    ## ASV5    767      5 OTUs
    ## ASV6    654      6 OTUs

On fait un graphique avec la fonction ggplot. On prend readsumsdf et on
met en x les attributions et en y le nombres de reads des séquences.
geom_bar est le nombre de séquences différentes. ggtitle permet
d’afficher un titre choisi. L’échelle en y est en log10.

``` r
ggplot(readsumsdf, aes(x = sorted, y = nreads)) +
  geom_bar(stat = "identity") +
  ggtitle("Total number of reads") +
  scale_y_log10() +
  facet_wrap(~type, nrow = 1, scales = "free")
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

On s’assure que l’effort d’échantillonage est le même. On définit un
minimum de lecture dans un échantillon.

``` r
set.seed(10000)
min(rowSums(physeq@otu_table@.Data))
```

    ## [1] 837

On crée une variable physeq_rar pour l’échantillonage soit de 800
lectures par échanillon.

``` r
physeq_rar <- rarefy_even_depth(physeq, sample.size = 800)
```

    ## You set `rngseed` to FALSE. Make sure you've set & recorded
    ##  the random seed of your session for reproducibility.
    ## See `?set.seed`

    ## ...

    ## 5OTUs were removed because they are no longer 
    ## present in any sample after random subsampling

    ## ...

``` r
rowSums(physeq_rar@otu_table@.Data)
```

    ## S11B  S1B  S2B  S2S  S3B  S3S  S4B  S4S  S5B  S5S  S6B  S6S  S7B  S7S  S8B  S8S 
    ##  800  800  800  800  800  800  800  800  800  800  800  800  800  800  800  800 
    ##  S9B  S9S 
    ##  800  800

``` r
physeq
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 213 taxa and 18 samples ]
    ## sample_data() Sample Data:       [ 18 samples by 21 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 213 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 213 tips and 212 internal nodes ]
    ## refseq()      DNAStringSet:      [ 213 reference sequences ]

``` r
physeq_rar
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 208 taxa and 18 samples ]
    ## sample_data() Sample Data:       [ 18 samples by 21 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 208 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 208 tips and 207 internal nodes ]
    ## refseq()      DNAStringSet:      [ 208 reference sequences ]

Deuxième méthode : nous utilisons des ratios. Pour cela, il faut
transformer nos données en ratios, ce qui nous permettra d’utiliser des
méthodes statistiques plus fidèles à la réalité. Pour cela, on crée un
objet tmp, qui prendra la fonction zCompositions::cmultRepl, avec method
= CZM (count zero multiplicative), label = 0 (par défaut), z.warning = 1
utilisé pour supprimer des colonnes ou des lignes incluant un excès de 0
ou de valeurs non observées. On effectue la transformation des données.
L’objet physeq_clr_asv prend la fonction apply, utilisée pour appliquer
la transformation à chaque ligne de la matrice tmp. La transformation
consiste à prendre le log de chaque valeur et à soustraire la moyenne du
log des valeurs de la ligne respective. Cela aide à centrer les données
et à les rendre appropriées pour certaines analyses statistiques.

``` r
tmp <- zCompositions::cmultRepl(physeq@otu_table,
                                method = "CZM",
                                label = 0,
                                z.warning = 1)
physeq_clr_asv <- apply(tmp, 1, function(x) log(x) - mean(log(x)))
```

On crée un nouvel objet physeq_clr en lui attribuant physeq. Ensuite, on
remplace la table OTU de cet objet par la matrice de données
transformées physeq_clr_asv après l’avoir transposée et convertie en une
matrice.

``` r
physeq_clr <- physeq
otu_table(physeq_clr) <- otu_table(t(physeq_clr_asv),
                                   taxa_are_rows = FALSE)
data.frame(physeq_clr@otu_table@.Data[1:5, 1:10])
```

    ##          ASV1       ASV2     ASV3       ASV4     ASV5     ASV6       ASV7
    ## S11B 5.172800  3.6295018 4.853277  4.6591212 4.876534 4.099505  4.4536772
    ## S1B  4.630559 -0.6264429 3.561361 -0.6264429 4.357692 4.297068 -0.6264429
    ## S2B  4.065517 -0.7557464 3.859665  3.0123670 4.041986 4.255561 -0.7557464
    ## S2S  5.042825  4.8740037 4.738829  2.8930022 5.003215 4.169296  3.9916145
    ## S3B  4.440233 -0.6954498 3.828432 -0.6954498 4.254516 4.653155 -0.6954498
    ##            ASV8       ASV9      ASV10
    ## S11B  3.9369865  4.1241980 -0.6524920
    ## S1B  -0.6264429  3.7217036  4.4863097
    ## S2B  -0.7557464 -0.7557464  4.0655169
    ## S2S   4.6847617  4.2367369 -0.6558837
    ## S3B  -0.6954498 -0.6954498  4.4057470

3.  On visualise l’abondance relative des organismes à des rangs
    taxonomiques spécifiques. On agglomère les données au niveau de la
    famille puis on transforme en abondance relative via la fonction
    transform_sample_counts. On filtre les taxons à faible abondance et
    on trie le bloc de données par ordre alphabétique par phylum.

``` r
physeq_phylum <- physeq_rar %>%
  tax_glom(taxrank = "Family") %>%                    
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() %>%                                       
  filter(Abundance > 0.02) %>%                        
  arrange(Family)                                      

head(physeq_phylum)
```

    ##    OTU Sample Abundance SampName   Geo Description groupe Pres PicoEuk Synec
    ## 1 ASV7    S4S 0.2483311      S4S North     North4S    NBS   78    2220  3130
    ## 2 ASV7    S4B 0.2129784      S4B North     North4B    NBF   78    2220  3130
    ## 3 ASV7    S5S 0.1833085      S5S North     North5S    NBS    0    1620 56555
    ## 4 ASV7    S2B 0.1776119      S2B North     North2B    NBF   59     890 25480
    ## 5 ASV7    S3S 0.1745283      S3S North     North3S    NBS    0     715 26725
    ## 6 ASV7    S5B 0.1741214      S5B North     North5B    NBF   42    1620 55780
    ##   Prochloro NanoEuk Crypto SiOH4   NO2   NO3   NH4   PO4    NT    PT   Chla
    ## 1     29835    2120    235 2.457 0.099 1.087 0.349 0.137 8.689 3.129 0.0000
    ## 2     29835    2120    235 2.457 0.099 1.087 0.349 0.137 8.689 3.129 0.0000
    ## 3     22835    2560    945 2.669 0.136 0.785 0.267 0.114 9.146 3.062 0.0000
    ## 4     16595     670    395 2.592 0.105 1.125 0.328 0.067 9.378 3.391 0.0000
    ## 5     16860     890    200 1.656 0.098 0.794 0.367 0.095 7.847 2.520 0.0000
    ## 6     23795    2555   1355 2.028 0.103 1.135 0.216 0.128 8.623 3.137 0.0102
    ##         T       S Sigma_t  Kingdom         Phylum               Class
    ## 1 18.8515 37.4542 26.9415 Bacteria Proteobacteria Alphaproteobacteria
    ## 2 18.8515 37.4542 26.9415 Bacteria Proteobacteria Alphaproteobacteria
    ## 3 24.1789 38.3213 26.1065 Bacteria Proteobacteria Alphaproteobacteria
    ## 4 22.6824 37.6627 26.0521 Bacteria Proteobacteria Alphaproteobacteria
    ## 5 22.5610 37.5960 26.0332 Bacteria Proteobacteria Alphaproteobacteria
    ## 6 24.1905 38.3192 26.1037 Bacteria Proteobacteria Alphaproteobacteria
    ##              Order                  Family
    ## 1 Rhodospirillales AEGEAN-169 marine group
    ## 2 Rhodospirillales AEGEAN-169 marine group
    ## 3 Rhodospirillales AEGEAN-169 marine group
    ## 4 Rhodospirillales AEGEAN-169 marine group
    ## 5 Rhodospirillales AEGEAN-169 marine group
    ## 6 Rhodospirillales AEGEAN-169 marine group

On visualise la répartition hiérarchique des phylums grâce à la fonction
treemap::treemap. fontside.labels permet de donner la taille des
étiquettes. fontcolor.labels permet de donner la couleur des étiquettes.
fontface.labels permet de définir une police aux étiquettes ainsi que
normal, gra, italique. align.labels permet de choisir l’alignement et la
position des étiquettes. overlap.labels permet de déterminer le
chevauchement des étiquettes (ici, la valeur par défaut 0,5 signifie que
les étiquettes de niveau inférieur sont imprimées si les autres
étiquettes ne se chevauchent pas sur plus de 0,5 fois leur taille de
zone). inflate.labels permet d’obtenir des étiquettes plus grandes si le
rectangle les contenant est grand. border.col permet de définir la
couleur des bordures de séparations.

``` r
treemap::treemap(physeq_phylum, index=c("Class", "Family"), vSize="Abundance", type="index",
        fontsize.labels=c(15,12),                
        fontcolor.labels=c("white","black"),    
        fontface.labels=c(2,1),                 
        align.labels=list(
          c("center", "center"), 
          c("left", "bottom")),                 
        overlap.labels=0.5,                     
        inflate.labels=F, 
        border.col=c("black","white"),          
        border.lwds=c(4,2),
        fontsize.title=12
)
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

On attribue à tmp la fonction transform_sample_counts qui transforme les
dénombrements d’échantillons dans une abondance taxonomique. Le
dénombrement pour chaque échantillon sera transformé individuellement.
La fonction ggplot permet de visualiser l’objet tmp précédemment
réalisé.

``` r
tmp <- transform_sample_counts(physeq,function(x) {x/sum(x)} ) %>%
  psmelt() %>%
  group_by(Family, Class) %>%
  summarise(abundance = sum(Abundance)) %>%
  na.omit()
```

    ## `summarise()` has grouped output by 'Family'. You can override using the
    ## `.groups` argument.

``` r
ggplot(tmp,aes(area=abundance,label=Family,fill=Class,subgroup=Class))+
  treemapify::geom_treemap()+
  treemapify::geom_treemap_subgroup_border() +
  treemapify::geom_treemap_subgroup_text(place = "centre",
                                         grow = T,
                                         alpha = 0.5,
                                         colour = "black",
                                         fontface = "italic",
                                         min.size = 0) +
  treemapify::geom_treemap_text(colour = "white",
                                place = "topleft",
                                reflow = TRUE)+
  theme(legend.position="none")
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

Nous pouvons ensuite sauvegarder notre figure dans le dossier
output_beta.

``` r
ggsave(here::here(output_beta,"treemap_treemapify.pdf"))
```

    ## Saving 7 x 5 in image

On reproduit les mêmes étapes avec notre table de composition en ASV
(physeq_phylum).

``` r
ggplot(physeq_phylum, aes(x = Sample, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity") +
  ylab("Relative Abundance (Family > 2%)") +
  scale_y_continuous(expand = c(0,0)) + 
  ggtitle("Community composition") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, size = 10,
                                   hjust = 0.5, vjust = 0.8),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())  
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
ggsave(here::here(output_beta, "asv_composition.pdf"))
```

    ## Saving 7 x 5 in image

4.  Ensuite, on calcule une matrice de distance grâce à l’indice de
    Jaccard. La distance de Jaccard est utilisée pour mesurer la
    similarité / dissimilarité de deux échantillons. Binary = TRUE
    signifie que les données seront lues comme absentes ou présentes
    (avec des 0 et des 1 dans la matrice respectivement). Cependant, en
    faisant cela, on risque des valeurs négatives donc on effectue la
    racine carrée de la matrice. Cela permet la préparation de la
    matrice de distance à la PCoA (analyse en coordonnées principales).

``` r
physeq_rar_jaccard <- phyloseq::distance(physeq_rar,
                                         method = "jaccard",
                                         binary = TRUE)
physeq_rar_jaccard <- sqrt(physeq_rar_jaccard)
```

Phylogénétique compositionnelle (unifrac non pondéré) : on vérifie si
l’arbre est bien enraciné et on calcule les distances Unifrac. unifracs
est une liste contenant 5 matrices de distance correspondant à l’UniFrac
pondéré (d_1), non pondéré (d_UW), ajusté à la variance (d_VAW), au
GUniFrac avec alpha = 0 et alpha = 0,5. GUniFrac permet de calculer la
matrice de distance unifracs. Elle permet de calculer les distances
UniFrac généralisées pour mesurer la dissemblance entre les communautés
microbiennes. Il prend en entrée le tableau OTU, l’arbre phylogénétique
et un vecteur de valeurs alpha (0, 0,5, 1) pour pondérer les distances
UniFrac.

``` r
ape::is.rooted(physeq_rar@phy_tree)
```

    ## [1] TRUE

``` r
unifracs <- GUniFrac::GUniFrac(physeq_rar@otu_table@.Data, physeq_rar@phy_tree, alpha=c(0, 0.5, 1))$unifracs
physeq_rar_du <- unifracs[, , "d_UW"]
```

L’objet tmp prend les valeurs normalisées de chaque échantillon. La
fonction transform_sample_counts prend chaque décompte d’un échantillon
et le divise par la somme des décomptes de cet échantillon. Cela nous
permet de calculer la dissimilarité de Bray-Curtis qui permet de mesurer
la dissimilarité entre des échantillons grâce à leur composition et leur
abondance en espèces. La fonction phyloseq::distance permet de calculer
la matrice de distance et (bray = méthode de calcul).

``` r
tmp <- transform_sample_counts(physeq,function(x) {x/sum(x)} )
physeq_rar_bray <- phyloseq::distance(tmp, method = "bray")
```

``` r
physeq_rar_dw <- unifracs[, , "d_1"] 
```

On peut calculer ces distances directement dans phyloseq en utilisant la
fonction dist.cal. Ici, on va faire une boucle qui parcourt chaque
méthode de distance.

``` r
dist_methods <- unlist(distanceMethodList)
data.frame(position = seq_along(dist_methods),
           dist_methods)
```

    ##             position dist_methods
    ## UniFrac1           1      unifrac
    ## UniFrac2           2     wunifrac
    ## DPCoA              3        dpcoa
    ## JSD                4          jsd
    ## vegdist1           5    manhattan
    ## vegdist2           6    euclidean
    ## vegdist3           7     canberra
    ## vegdist4           8         bray
    ## vegdist5           9   kulczynski
    ## vegdist6          10      jaccard
    ## vegdist7          11        gower
    ## vegdist8          12     altGower
    ## vegdist9          13     morisita
    ## vegdist10         14         horn
    ## vegdist11         15    mountford
    ## vegdist12         16         raup
    ## vegdist13         17     binomial
    ## vegdist14         18         chao
    ## vegdist15         19          cao
    ## betadiver1        20            w
    ## betadiver2        21           -1
    ## betadiver3        22            c
    ## betadiver4        23           wb
    ## betadiver5        24            r
    ## betadiver6        25            I
    ## betadiver7        26            e
    ## betadiver8        27            t
    ## betadiver9        28           me
    ## betadiver10       29            j
    ## betadiver11       30          sor
    ## betadiver12       31            m
    ## betadiver13       32           -2
    ## betadiver14       33           co
    ## betadiver15       34           cc
    ## betadiver16       35            g
    ## betadiver17       36           -3
    ## betadiver18       37            l
    ## betadiver19       38           19
    ## betadiver20       39           hk
    ## betadiver21       40          rlb
    ## betadiver22       41          sim
    ## betadiver23       42           gl
    ## betadiver24       43            z
    ## dist1             44      maximum
    ## dist2             45       binary
    ## dist3             46    minkowski
    ## designdist        47          ANY

``` r
dist_methods <- dist_methods[c(1, 2, 10, 8)]
dist_methods
```

    ##   UniFrac1   UniFrac2   vegdist6   vegdist4 
    ##  "unifrac" "wunifrac"  "jaccard"     "bray"

On calcule la matrice de distance en utilisant la méthode de distance
actuelle, qu’on attribue ensuite à l’objet iDist. On fait l’ordination
PCoA grâce à cette matrice, et on l’attribue à iMDS. À partir de cela,
on peut tracer l’ordination qu’on attribue à p puis on la colorie via
l’attribu Geo. On n’oublie pas le titre du tracé qui indique la méthode
de distance utilisée. On enregistre chaque tracé dans une liste “plist”
qui sera nommé par le nom de la méthode de distance utilisée.

``` r
plist <- vector("list")

for(i in dist_methods){
  iDist <- phyloseq::distance(physeq_rar, method = i)
  iMDS <- ordinate(physeq_rar, "MDS", distance = iDist)
  p <- NULL
  p <- plot_ordination(physeq_rar, iMDS, color= "Geo")
  p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
  plist[[i]] = p 
}
```

On combine les résulats…

``` r
df <- plyr::ldply(plist, function(x) x$data)
head(df)
```

    ##       .id      Axis.1      Axis.2 SampName   Geo Description groupe Pres
    ## 1 unifrac  0.09023445  0.06150644     S11B South     South5B    SGF   35
    ## 2 unifrac -0.21048836 -0.19946687      S1B North     North1B    NBF   52
    ## 3 unifrac -0.21001002 -0.08655455      S2B North     North2B    NBF   59
    ## 4 unifrac  0.12583068  0.07022248      S2S North     North2S    NBS    0
    ## 5 unifrac -0.31465014 -0.06077941      S3B North     North3B    NBF   74
    ## 6 unifrac -0.16616937  0.01827175      S3S North     North3S    NBS    0
    ##   PicoEuk Synec Prochloro NanoEuk Crypto SiOH4   NO2   NO3   NH4   PO4    NT
    ## 1    5370 46830       580    6010   1690 3.324 0.083 0.756 0.467 0.115 9.539
    ## 2     660 32195     10675     955    115 1.813 0.256 0.889 0.324 0.132 9.946
    ## 3     890 25480     16595     670    395 2.592 0.105 1.125 0.328 0.067 9.378
    ## 4     890 25480     16595     670    395 3.381 0.231 0.706 0.450 0.109 8.817
    ## 5     835 13340     25115    1115    165 1.438 0.057 1.159 0.369 0.174 8.989
    ## 6     715 26725     16860     890    200 1.656 0.098 0.794 0.367 0.095 7.847
    ##      PT   Chla       T       S Sigma_t
    ## 1 4.138 0.0182 23.0308 38.9967 26.9631
    ## 2 3.565 0.0000 22.7338 37.6204 26.0046
    ## 3 3.391 0.0000 22.6824 37.6627 26.0521
    ## 4 3.345 0.0000 22.6854 37.6176 26.0137
    ## 5 2.568 0.0000 21.5296 37.5549 26.2987
    ## 6 2.520 0.0000 22.5610 37.5960 26.0332

…puis on change le bloc de données pour donner à la première colonne le
nom de distance. On utilise df comme source de donnés, ce qui nous
permet de faire le plot. La variable aes permet de spécifier
l’esthétique du plot, avec toujours GEO en couleur. On peut ajouter des
points au tracé, d’une taille de 3 et d’une transparence de 0.5 qui
permet de créer un nuage de point. On va utiliser le thème noir et blanc
(theme_bw) puis créer des graphiques de plus petites tailles basés sur
la variable distance, sachant que chaque graphique aura sa propre
échelle (scales=free). On n’oublie pas le titre du plot.

``` r
names(df)[1] <- "distance"

ggplot(df, aes(Axis.1, Axis.2, color = Geo)) +
  geom_point(size=3, alpha=0.5) +
  theme_bw() +
  facet_wrap(~distance, scales="free") +
  ggtitle("PCoA (MDS) on various distance metrics")
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

5.  On peut maintenant faire le regroupement hiérarchique. On va se
    baser sur la distance d’Aitchison. D’abord, on va créer un objet qui
    va accueillir le calcul de la matrice de distance, avec la méthode
    euclidienne.

``` r
physeq_clr_dist <- phyloseq::distance(physeq_clr, method = "euclidean")
```

On veut effectuer un regroupement hiérarchique ascendant sur la matrice
de distance. hclust permet de calculer une structure hiérarchique de
clusters sous forme de dendrogrammes (arbres binaires où les feuilles
représentent les observations individuelles et les nœuds internes
représentent les groupes ou clusters). Ici, chaque ligne représente une
méthode particulière : simple agrégation, agrégation complète,
arthimétique ou par parties. par(mfrow=) permet d’obtenir 2x2
graphiques.

``` r
spe_single <- hclust(physeq_clr_dist, method = "single")
spe_complete <- hclust(physeq_clr_dist, method = "complete")
spe_upgma <- hclust(physeq_clr_dist, method = "average")
spe_ward <- hclust(physeq_clr_dist, method = "ward.D")

par(mfrow = c(2, 2))
plot(spe_single, main = "single")
plot(spe_complete, main = "complete")
plot(spe_upgma, main = "UPGMA")
plot(spe_ward, main = "ward")
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

Ensuite, on fait une matrice cophénétique (matrice représentant les
distances cophénétiques entre toutes les paires d’objets). Une
corrélation r de Pearson, appelée corrélation cophénétique dans ce
contexte, peut être calculée entre la matrice de dissimilarité originale
et la matrice cophénétique. La méthode présentant la corrélation
cophénétique la plus élevée peut être considérée comme celle qui a
produit le meilleur modèle de regroupement pour la matrice de distance.
On calcule la matrice cophénétique et la corrélation des quatre
résultats de clustering présentés ci-dessus, au moyen de la fonction
cophenetic() du package stats.

``` r
spe_single_coph <- cophenetic(spe_single)
cor(physeq_clr_dist, spe_single_coph)
```

    ## [1] 0.9447202

``` r
spe_complete_coph <- cophenetic(spe_complete)
cor(physeq_clr_dist, spe_complete_coph)
```

    ## [1] 0.8609329

``` r
spe_upgma_coph <- cophenetic(spe_upgma)
cor(physeq_clr_dist, spe_upgma_coph)
```

    ## [1] 0.958006

``` r
spe_ward_coph <- cophenetic(spe_ward)
cor(physeq_clr_dist, spe_ward_coph)
```

    ## [1] 0.9044309

Pour illustrer la relation entre une matrice de distance et un ensemble
de matrices cophénétiques, on peut traçer les distances originales par
rapport aux distances cophénétiques.

``` r
plot_coph_cor <- function(cophenetic_distance, hclust_type){
  cor_res <- round(cor(physeq_clr_dist, cophenetic_distance),3)
  plot(x = physeq_clr_dist,
     y = cophenetic_distance,
     xlab = "Aitchison distance",
     ylab = "Cophenetic distance",
     xlim = c(10, 35), ylim = c(10, 35),
     main = c(hclust_type, paste("Cophenetic correlation ", cor_res)))
  abline(0, 1)
}

par(mfrow=c(2,2))

plot_coph_cor(cophenetic_distance = spe_complete_coph,
              hclust_type = "Single linkage")

plot_coph_cor(cophenetic_distance = spe_complete_coph,
              hclust_type = "Complete linkage")

plot_coph_cor(cophenetic_distance = spe_upgma_coph,
              hclust_type = "Average linkage")

plot_coph_cor(cophenetic_distance = spe_ward_coph,
              hclust_type = "Ward linkage")
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

Pour interpréter et comparer les résultats du clustering, il faut des
clusters interprétables. Les valeurs du niveau de fusion d’un
dendrogramme sont les valeurs de dissimilarité où se produit une fusion
entre deux branches d’un dendrogramme. Tracer les valeurs du niveau de
fusion peut aider à définir les niveaux de coupe.

``` r
par(mfrow = c(1, 1))

plot(x = spe_upgma$height,
     y = phyloseq::nsamples(physeq_clr):2,
     type = "S",
     main = "Fusion levels - Aitchison - Average",
     ylab = "k (number of cluster)",
     xlab = "h (node height)")

text(x = spe_upgma$height,
     y = phyloseq::nsamples(physeq_clr):2,
     labels = phyloseq::nsamples(physeq_clr):2,
     col = "red",
     cex = 0.8)
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

On installe le package NbClust qui calculera, avec un seul appel de
fonction, 24 indices pour confirmer le bon nombre de clusters dans
l’ensemble de données. NbClust confirme l’identification de deux groupes
d’échantillons. Nous revenons sur le dendrogramme et le découpons aux
distances correspondantes.

``` r
install.packages("NbClust", lib = ".")
library("NbClust", lib.loc = ".")
nclust <- nb_clust_all(data = t(physeq_clr_asv), seed = 1000)
```

    ## [1] "Trying kl index..."
    ## [1] "Trying ch index..."
    ## [1] "Trying hartigan index..."
    ## [1] "Trying scott index..."
    ## [1] "Trying cindex index..."
    ## [1] "Trying db index..."
    ## [1] "Trying silhouette index..."
    ## [1] "Trying duda index..."
    ## [1] "Trying pseudot2 index..."
    ## [1] "Trying beale index..."
    ## [1] "Trying ratkowsky index..."
    ## [1] "Trying ball index..."
    ## [1] "Trying ptbiserial index..."
    ## [1] "Trying gap index..."
    ## [1] "Trying frey index..."
    ## [1] "Trying mcclain index..."
    ## [1] "Trying gamma index..."
    ## [1] "Trying gplus index..."
    ## [1] "Trying tau index..."
    ## [1] "Trying dunn index..."
    ## [1] "Trying hubert index..."

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

    ## *** : The Hubert index is a graphical method of determining the number of clusters.
    ##                 In the plot of Hubert index, we seek a significant knee that corresponds to a 
    ##                 significant increase of the value of the measure i.e the significant peak in Hubert
    ##                 index second differences plot. 
    ##  
    ## [1] "Trying sdindex index..."
    ## [1] "Trying dindex index..."

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-32-2.png)<!-- -->

    ## *** : The D index is a graphical method of determining the number of clusters. 
    ##                 In the plot of D index, we seek a significant knee (the significant peak in Dindex
    ##                 second differences plot) that corresponds to a significant increase of the value of
    ##                 the measure. 
    ##  
    ## [1] "Trying sdbw index..."
    ## Based on a number of criteria, we will select 2 clusters.

On découpe le dendrogramme créé par clustering hiérarchique (méthode de
liaison UPGMA) en groupes k, puis on examine la composition de ces
groupes à l’aide de la fonction table.

``` r
k <- 2
spe_upgma_clust <- cutree(tree = spe_upgma, k = k)
table(spe_upgma_clust)
```

    ## spe_upgma_clust
    ##  1  2 
    ## 12  6

``` r
spe_upgma_clust2 <- data.frame(UPGMA_clusters = spe_upgma_clust)
```

On trace ensuite le dendrogramme avec les groupes.

``` r
plot(spe_upgma,
     hang = -1,
     ylab = "Height",
     main="Aitchison distance - UPGMA")

rect.hclust(spe_upgma,
            k = k,
            border = 2:6,
            cluster = spe_upgma_clust)

legend("topright",
       paste("Cluster", 1:k),
       pch = 22,
       col = 2:(k + 1),
       bty = "n")
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

Il existe plusieurs façons de mesurer la robustesse d’un algorithme de
clustering. Trois mesures couramment utilisées sont l’indice Dunn
(rapport de la plus petite distance inter-cluster à la plus grande
distance intra-cluster), l’indice Davis-Bouldin et l’indice Silhoutte.
Un DI élevé signifie un meilleur regroupement car les observations de
chaque cluster sont plus rapprochées. on utilise cluster.stats() pour
calculer l’index Dunn.

``` r
cs <- fpc::cluster.stats(d = physeq_clr_dist,
                         clustering = spe_upgma_clust)

cs$dunn
```

    ## [1] 0.9231545

L’indice de Dunn est élevé, ce qui indique un bon regroupement des
échantillons. Maintenant que nous avons identifié deux groupes
d’échantillons en fonction de la composition de leur communauté
microbienne, nous souhaiterons peut-être examiner quels clades
microbiens ou ASV sont enrichis dans chacun des groupes.

On va combiner le Z-score et le clustering. Pour cela, on normalise les
valeurs (on centre autour de la moyenne et on réduit). C’est la
comparaison d’une valeur observée d’un échantillon à la moyenne de la
population. Elle répond donc à la question , à quelle distance de la
moyenne de la population se trouve un score pour un échantillon donné.
Les scores sont donnés en SD par rapport à la moyenne de la population.

On transforme les comptages de chaque échantillon en pourcentage
(transform_sample_counts). Ensuite, on calcule la somme des comptages
pour chaque taxon (taxa_sums). On les trie dans l’ordre décroissant et
on ne prend que les 30 plus abondants. On crée un objet selection30 qui
prend les 30 taxons les plus abondants.

``` r
pourcentS <- phyloseq::transform_sample_counts(physeq_rar, function(x) x/sum(x) * 100)
mytop30 <- names(sort(phyloseq::taxa_sums(pourcentS), TRUE)[1:30])
selection30 <- phyloseq::prune_taxa(mytop30, pourcentS)
selection30
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 30 taxa and 18 samples ]
    ## sample_data() Sample Data:       [ 18 samples by 21 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 30 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 30 tips and 29 internal nodes ]
    ## refseq()      DNAStringSet:      [ 30 reference sequences ]

On extrait la table des ASV de l’objet selection30, puis on extrait les
données associées qu’on stocke dans l’objet selection30_sample.

``` r
selection30_asv <- phyloseq::otu_table(selection30)
selection30_sample <- phyloseq::sample_data(selection30)
rownames(selection30_asv)
```

    ##  [1] "S11B" "S1B"  "S2B"  "S2S"  "S3B"  "S3S"  "S4B"  "S4S"  "S5B"  "S5S" 
    ## [11] "S6B"  "S6S"  "S7B"  "S7S"  "S8B"  "S8S"  "S9B"  "S9S"

On crée un nouvel objet qui prend comme valeur la concaténation des
colonnes SampName et Description avec \_ comme séparateur. Ensuite,
l’objet heat prend comme valeur la transposition de la table ASV
(échange des lignes par les colonnes) qui a centré puis réduit les
valeurs. La matrice heat est ensuite convertie en data.frame. On affiche
les 6 premières lignes avec la fonction head.

``` r
sample_new_names <- paste(selection30_sample$SampName,
                          selection30_sample$Description,
                          sep = "_")

heat <- t(base::scale(selection30_asv))
head(data.frame(heat))
```

    ##            S11B         S1B        S2B        S2S        S3B         S3S
    ## ASV1  1.0670101 -0.36085474 -0.8368097  0.5070631 -0.6688256  0.08710287
    ## ASV2 -0.3822681 -0.72549212 -0.7254921  0.5166521 -0.7254921 -0.59474010
    ## ASV3  1.4657223 -1.12279860 -0.5254476  0.6692543 -0.3761099 -2.16816282
    ## ASV4  0.2776466 -1.18129127 -0.9019202 -0.8087965 -1.1812913 -0.77775527
    ## ASV5  1.1642633  0.31674810  0.2397013  1.3954038  0.3552715  0.20117785
    ## ASV6  0.4514863 -0.01289961  0.7417276  0.3934381  1.4963548 -1.81239520
    ##             S4B        S4S         S5B        S5S         S6B         S6S
    ## ASV1 -1.7327249 -0.3608547  1.48697039  2.2149015  1.48697039 -0.47284414
    ## ASV2 -0.6110841 -0.6764601 -0.49667608 -0.7254921  0.38590007  3.34416457
    ## ASV3 -0.6747854 -0.7245646  1.06748833 -1.0232401 -0.02765514  0.02212411
    ## ASV4 -0.3121368 -1.1812913  1.67450193  1.1778422 -0.68463158 -0.56046666
    ## ASV5 -0.3766734  0.5864120 -1.30123544  0.2782247  1.43392721 -1.30123544
    ## ASV6  0.7997758 -1.8123952  0.04514863 -1.8123952 -0.59338206 -0.59338206
    ##             S7B         S7S        S8B        S8S        S9B          S9S
    ## ASV1 -0.8928044  0.03110817 -0.6128309 -0.3888521 -0.5568362  0.003110817
    ## ASV2  0.9742842 -0.18614003  0.4349321 -0.5293641  0.4022441  0.320524054
    ## ASV3  0.6194751 -0.22677213  0.5696958  1.8639563 -0.1769929  0.768812836
    ## ASV4  0.8363887  1.02263609  1.2088835  1.1778422  0.8053475 -0.591507891
    ## ASV5 -1.3012354 -1.30123544 -1.3012354  0.3552715  1.2798335 -0.723384173
    ## ASV6 -0.1870443  0.10319688  0.7417276  0.2192934  0.5095346  1.322210022

On visualise les données sous la forme de carte thermique (fonction
Heatmap). Ainsi, on visualise la matrice heat, on ajuste la taille de la
police des noms d’échantillons puis on désactive le clustering pour
garder l’ordre original des colonnes. On personnalise la légende de la
carte thermique : verticale, avec Z-scores en titre, la largeur et la
hauteur de 0.5 cm et 3 cm respectivement.

``` r
ComplexHeatmap::Heatmap(
  heat,
  row_names_gp = grid::gpar(fontsize = 6),
  cluster_columns = FALSE,
  heatmap_legend_param = list(direction = "vertical",
                              title = "Z-scores", 
                              grid_width = unit(0.5, "cm"),
                              legend_height = unit(3, "cm"))
)
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

On extrait la table de selection30 à l’aide de tax_table. On stocke les
variables dans l’objet taxon. On change le nom des colonnes en combinant
les noms de lignes de taxon, les infos du phylum et les infos de la
famille, qui sont séparés par des \_. Puis on remplace les noms de
colonnes actuels de la table selection30_asv par les nouveaux noms de
myname.

``` r
taxon <- phyloseq::tax_table(selection30) |>
  as.data.frame()
myname <- paste(rownames(taxon), taxon$Phylum, taxon$Family, sep="_")
colnames(selection30_asv) <- myname
```

On prend la matrice selection30_asv qu’on centre-réduit et qu’on
transpose dans heat. On crée un objet my_top_annotation qui ajoute des
infos supplémentaires au-dessus de la carte thermique. On la visualise
avec des paramètres personnalisés (taille de police,…).

``` r
heat <- t(scale(selection30_asv))

my_top_annotation <- ComplexHeatmap::anno_block(gp = grid::gpar(fill =c(3,4)),
                                               labels = c(1, 2),
                                               labels_gp = grid::gpar(col = "white",
                                                                      fontsize = 10))

ComplexHeatmap::Heatmap(
  heat,
  row_names_gp = grid::gpar(fontsize = 6),
  cluster_columns =TRUE,
  heatmap_legend_param = list(direction = "vertical",
   title ="Z-scores",
   grid_width = unit(0.5, "cm"),
   legend_height = unit(4, "cm")),
  top_annotation = ComplexHeatmap::HeatmapAnnotation(foo = my_top_annotation),
  column_km = 2,
  column_names_gp= grid::gpar(fontsize = 6)
  )
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-41-1.png)<!-- -->

Un boxplot est généré pour chaque ligne de selection30_asv montrant la
distribution des abondances. Ces boxplots sont affichés comme une
annotation à gauche de la carte thermique principale
(my_box_plot_left_anno). De la même manière, on ajoute des informations
supplémentaires en haut (my_top_anno). Ensuite, on visualise les donnés
sous la forme de carte thermique avec toutes les informations
précedemment demandées.

``` r
boxplot <- ComplexHeatmap::anno_boxplot(t(selection30_asv), 
                                        which = "row",
                                        gp = grid::gpar(fill = "turquoise3"))

my_boxplot_left_anno <- ComplexHeatmap::HeatmapAnnotation(Abund = boxplot,
                                                          which = "row",
                                                          width = unit(3, "cm"))

my_top_anno <- ComplexHeatmap::anno_block(gp = grid::gpar(fill = c(3, 6)),
                                          labels = c("South", "North"),
                                          labels_gp = grid::gpar(col = "white",
                                                                fontsize = 10))

my_top_anno <- ComplexHeatmap::HeatmapAnnotation(foo = my_top_anno)

ComplexHeatmap::Heatmap(
  heat,
  row_names_gp = grid::gpar(fontsize = 7),
  left_annotation = my_boxplot_left_anno, 
  heatmap_legend_param = list(direction = "vertical",
                              title ="Z-scores",
                              grid_width = unit(0.5, "cm"),
                              legend_height = unit(3, "cm")),
  top_annotation = my_top_anno,
  column_km = 2,
  cluster_columns = TRUE,
  column_dend_side = "bottom",
  column_names_gp = grid::gpar(fontsize = 7)
  )
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-42-1.png)<!-- -->

6.  Analyse en composante principale (PCoA) : permet d’avoir un aperçu
    des relations entre les objets et les variables. On utilise là aussi
    la distance de Aitchinson. On commence par extraire la table de
    physeq_clr pour la transformer en data frame. L’objet ASVname prend
    la concaténation des noms de lignes actuels + famille et genre pour
    chaque taxon. Cela permet de renommer les lignes de physeq_clr_asv.
    La PCoA est réalisée sur physeq_clr_asv, les données associées à
    chaque échantillon de physeq_clr sont intégrées. On visualise les
    résulats par un screeplot, ce qui permet de déterminer le nombre de
    composantes à considérer.

``` r
tax_CLR <-  as.data.frame(tax_table(physeq_clr))
ASVname <- paste(rownames(tax_CLR), tax_CLR$Family, tax_CLR$Genus,sep="_")
rownames(physeq_clr_asv) <- ASVname
p <- PCAtools::pca(physeq_clr_asv,
                   metadata = data.frame(sample_data(physeq_clr)))
PCAtools::screeplot(p, axisLabSize = 18, titleLabSize = 22)
```

    ## Warning: Removed 2 rows containing missing values (`geom_line()`).

    ## Warning: Removed 2 rows containing missing values (`geom_point()`).

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-43-1.png)<!-- -->

On détermine le nombre optimal de composantes principales à retenir. On
fait un test de Horn (parallelPCA) qui détermine le nombre de
composantes à conversver en comparant la variance de chaque composante
principlae à celle obtenue à partir de jeux de donnés aléatoires. On
extrait ce nombre avec horn\$n. On fait un test de Elbow qui identifie
le point où l’ajout de composantes supplémentaires n’est pas
significatif (on atteint donc la stabilisation de la variance). On
affiche ensuite cet indice.

``` r
horn <- PCAtools::parallelPCA(physeq_clr_asv)
```

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

    ## Warning in check_numbers(x, k = k, nu = nu, nv = nv): more singular
    ## values/vectors requested than available

``` r
horn$n
```

    ## [1] 2

``` r
elbow <- PCAtools::findElbowPoint(p$variance)
elbow
```

    ## PC3 
    ##   3

On visualise les résultats obtenus. La fonction biplot affcihe les
scores des échantillons et les vecteurs de chargement des variables. Les
points représentant les échantillons seront étiquettés avec les noms des
échantillons contenus dans p\$meta…. Encore une fois, on applique le
code couleur Geo. La taille des points est de 5, on trace des axes de
référence avec hline et vline. La légende est plaxée à droite du biplot.

``` r
PCAtools::biplot(
  p,
  lab = p$metadata$SampName,
  colby = "Geo",
  pointSize = 5,
  hline = 0, vline = 0,
  legendPosition = "right"
)
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-45-1.png)<!-- -->

On fait la même chose mais en ajoutant des paramètres supplémentaires.
On ajoute les vecteurs de chargements (donnent une indication de
l’importance de la direction de chaque variable dans l’espace des CP).
La longueur des flèches est augmentée de 50%. La taille des étiquettes
pour le noms des vecteurs de chargement est de 3 et on les colore en
rouge. On n’affiche que les 3 vecteurs de chargements les plus
importants. Puis, on étiquette les échantillons grâce à leur valeur de
X.SampleID. On utilise encore la coloration Geo, des lignes vline et
hline et la légende à droite.

``` r
PCAtools::biplot(
  p, 
  showLoadings = TRUE,
  lengthLoadingsArrowsFactor = 1.5,
  sizeLoadingsNames = 3,
  colLoadingsNames = 'red4',
  ntopLoadings = 3,
  lab = p$metadata$X.SampleID,
  colby = "Geo",
  hline = 0, vline = 0,
  legendPosition = "right"
)
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-46-1.png)<!-- -->

On visualise les corrélations entre les principales composantes d’une
PCoA et certains variables environnementales ou donnés. components :
composantes définies par les résultats du test de Horn. metavars : liste
des variables environnementales pour lesquelles les corrélations seront
calculées. On utilise plein de couleurs. On ajuste la taille et le style
de la police des valeurs de corrélation ainsi que la position des
étiquettes des axes. On utilise la corrélation de Spearman. On utilise
des symboles pour indiquer la significativité des corrélations.

``` r
PCAtools::eigencorplot(
  p,
  components = PCAtools::getComponents(p, 1:horn$n),
  metavars = c('SiOH4','NO2','NO3','NH4','PO4',
              'NT','PT','Chla',"T", "S", "Sigma_t"),
  col = c('white', 'cornsilk1', 'gold',
          'forestgreen', 'darkgreen'),
  cexCorval = 1.2,
  fontCorval = 2,
  posLab = "all",
  rotLabX = 45,
  scale = TRUE,
  main = bquote(PC ~ Spearman ~ r^2 ~ environmental ~ correlates),
  plotRsquared = TRUE,
  corFUN = "spearman",
  corUSE = "pairwise.complete.obs",
  corMultipleTestCorrection = 'BH',
  signifSymbols = c("****", "***", "**", "*", ""),
  signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1)
)
```

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

    ## Warning in cor.test.default(xvals[, i], yvals[, j], use = corUSE, method =
    ## corFUN): Cannot compute exact p-value with ties

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-47-1.png)<!-- -->

L’analyse PCoA est réalisée sur les distances de Bray-Curtis à partir de
physeq_rar_bray. Les coordonnées des deux premiers axes sont extraites
et stockées dans pcoa_coord. Le data frame hull contient les coordonnées
PCoA et les données associées. On donne deux couleurs pour les
catégories North et South. On crée les enveloppes convexes avec le
package dplyr. On affiche les premières lignes de hull_data pour
vérifier sa structure.

``` r
pcoa_asv <- ape::pcoa(physeq_rar_bray)
pcoa_coord <- pcoa_asv$vectors[, 1:2]
hull <- data.frame("Axis.1" = pcoa_coord[, 1],
                   "Axis.2" = pcoa_coord[, 2],
                   "sample" = as.data.frame(sample_data(physeq_rar@sam_data)))

hull_col <- c("#a65628","#1919ff")
names(hull_col) <- c("North","South")

hull_data <- hull %>%
  dplyr::group_by(sample.Geo) %>%
  dplyr::slice(chull(Axis.1,Axis.2)) %>%
  dplyr::mutate(color = hull_col[sample.Geo])

head(hull_data)
```

    ## # A tibble: 6 × 24
    ## # Groups:   sample.Geo [1]
    ##     Axis.1  Axis.2 sample.SampName sample.Geo sample.Description sample.groupe
    ##      <dbl>   <dbl> <chr>           <chr>      <chr>              <chr>        
    ## 1  0.242   -0.0992 S2S             North      North2S            NBS          
    ## 2 -0.403   -0.130  S2B             North      North2B            NBF          
    ## 3 -0.455   -0.0922 S3B             North      North3B            NBF          
    ## 4 -0.471    0.0176 S3S             North      North3S            NBS          
    ## 5  0.00454  0.407  S5S             North      North5S            NBS          
    ## 6  0.102    0.327  S5B             North      North5B            NBF          
    ## # ℹ 18 more variables: sample.Pres <int>, sample.PicoEuk <int>,
    ## #   sample.Synec <int>, sample.Prochloro <int>, sample.NanoEuk <int>,
    ## #   sample.Crypto <int>, sample.SiOH4 <dbl>, sample.NO2 <dbl>,
    ## #   sample.NO3 <dbl>, sample.NH4 <dbl>, sample.PO4 <dbl>, sample.NT <dbl>,
    ## #   sample.PT <dbl>, sample.Chla <dbl>, sample.T <dbl>, sample.S <dbl>,
    ## #   sample.Sigma_t <dbl>, color <chr>

On visualise les résultats de la PCoA sous forme de graphique. Il montre
comment les échantillons se répartissent dans l’espace bidimensionnel
des deux premières composantes principales et comment ils se regroupent
en fonction des catégories “North” et “South”. Les enveloppes convexes
aident à visualiser la distribution globale de chaque catégorie dans cet
espace.

``` r
ggplot(data = hull, aes(x = Axis.1, y = Axis.2)) +
  geom_hline(yintercept = 0, colour = "lightgrey", linetype = 2) +
  geom_vline(xintercept = 0, colour = "lightgrey", linetype = 2) +
  geom_polygon(data = hull_data,
               aes(group = sample.Geo,
                   fill = sample.Geo),
               alpha = 0.3) + # add the convex hulls)
  scale_fill_manual(values = c("Darkgrey", "#1919ff")) +
  geom_point(data = hull,
             aes(color = sample.Geo,
                 size = sample.S),
             alpha = 0.7) +
  scale_color_manual(values = c("Darkgrey", "#1919ff")) +
  xlab(paste("PCo1 (", round(pcoa_asv$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa_asv$values$Relative_eig[2]*100, 1), "%)")) +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size = 14), # remove x-axis labels
        axis.title.y = element_text(size = 14), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-49-1.png)<!-- -->

On réalise une analyse NMDS (représente fidèlement les distances entre
les objets dans un espace à dimension réduite) sur les distances de
physeq_clr_dist en utilisant le package vegan. k=2 : deux dimensions.
trymax=100 on tente 100 configurations pour trouver la meilleure, puis
les résultats sont stockés dans physeq_clr_nmds. La fonction stressplot
permet de visualiser le niveau de stress de l’analyse, avec le niveau de
stress associé à cette configuration.

``` r
physeq_clr_nmds <- vegan::metaMDS(physeq_clr_dist, k=2, trymax=100)
```

    ## Run 0 stress 0.06963252 
    ## Run 1 stress 0.06963252 
    ## ... New best solution
    ## ... Procrustes: rmse 3.519632e-06  max resid 6.819366e-06 
    ## ... Similar to previous best
    ## Run 2 stress 0.08705507 
    ## Run 3 stress 0.06963252 
    ## ... Procrustes: rmse 1.984082e-06  max resid 3.686583e-06 
    ## ... Similar to previous best
    ## Run 4 stress 0.08747761 
    ## Run 5 stress 0.08909908 
    ## Run 6 stress 0.06963252 
    ## ... Procrustes: rmse 4.436839e-06  max resid 9.114362e-06 
    ## ... Similar to previous best
    ## Run 7 stress 0.06963252 
    ## ... Procrustes: rmse 2.967907e-06  max resid 6.047389e-06 
    ## ... Similar to previous best
    ## Run 8 stress 0.08747757 
    ## Run 9 stress 0.06963252 
    ## ... Procrustes: rmse 9.883376e-07  max resid 2.056696e-06 
    ## ... Similar to previous best
    ## Run 10 stress 0.06963252 
    ## ... Procrustes: rmse 5.333631e-07  max resid 1.187664e-06 
    ## ... Similar to previous best
    ## Run 11 stress 0.06963252 
    ## ... Procrustes: rmse 2.401433e-06  max resid 4.837389e-06 
    ## ... Similar to previous best
    ## Run 12 stress 0.08747757 
    ## Run 13 stress 0.06963252 
    ## ... Procrustes: rmse 1.821518e-06  max resid 3.14952e-06 
    ## ... Similar to previous best
    ## Run 14 stress 0.06963252 
    ## ... Procrustes: rmse 3.91536e-06  max resid 7.684323e-06 
    ## ... Similar to previous best
    ## Run 15 stress 0.08909787 
    ## Run 16 stress 0.06963252 
    ## ... Procrustes: rmse 3.690599e-06  max resid 7.475907e-06 
    ## ... Similar to previous best
    ## Run 17 stress 0.06963252 
    ## ... Procrustes: rmse 2.399131e-06  max resid 4.89376e-06 
    ## ... Similar to previous best
    ## Run 18 stress 0.06963252 
    ## ... Procrustes: rmse 7.790933e-06  max resid 1.562559e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.06963252 
    ## ... Procrustes: rmse 1.300643e-06  max resid 2.74642e-06 
    ## ... Similar to previous best
    ## Run 20 stress 0.08705507 
    ## *** Best solution repeated 13 times

``` r
vegan::stressplot(physeq_clr_nmds)
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-50-1.png)<!-- -->

On visualise les résultats de l’analyse NMDS pour un ensemble de
données. Le graphique affiche les points représentant les échantillons
dans un espace bidimensionnel défini par les deux axes principaux de
l’NMDS. De plus, il trace des enveloppes convexes (hulls) autour des
points associés à des catégories spécifiques (dans ce cas, “North” et
“South”), aidant à visualiser la distribution et la séparation de ces
catégories dans l’espace NMDS. Le niveau de stress de l’analyse NMDS,
qui évalue la qualité de la représentation bidimensionnelle, est
également affiché sur le graphique.

``` r
nmds_coord <- data.frame(physeq_clr_nmds$points)

hull <- data.frame("Axis.1" = nmds_coord[,1],
                   "Axis.2" = nmds_coord[,2],
                   "sample" = as.data.frame(sample_data(physeq_clr@sam_data)))

hull_col <- c("#a65628","#1919ff")
names(hull_col) <- c("North","South")

hull_data <- hull %>%
  dplyr::group_by(sample.Geo) %>%
  dplyr::slice(chull(Axis.1,Axis.2)) %>%
  dplyr::mutate(color = hull_col[sample.Geo])

ggplot(hull,aes(x = Axis.1, y = Axis.2)) +
  geom_hline(yintercept = 0, colour = "lightgrey", linetype = 2) + 
  geom_vline(xintercept = 0, colour = "lightgrey", linetype = 2) +
  geom_polygon(data = hull_data,
               aes(group = sample.Geo,
                   fill = sample.Geo),
               alpha = 0.3) + 
  scale_fill_manual(values = c("Darkgrey", "#1919ff")) +
  geom_point(data = hull,
             aes(color = sample.Geo,
                 size = sample.S),
             alpha = 0.7) +
  scale_color_manual(values = c("Darkgrey", "#1919ff")) +
  geom_text(data = hull_data,
            x = -0, y = -9,
            label = paste("Stress =", round(physeq_clr_nmds$stress, 2)),
            colour = "Black",
            size = 5)  +
  xlab(paste("MDS1")) +
  ylab(paste("MDS2")) +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-51-1.png)<!-- -->

On évalue la correspondance entre les coordonnées NMDS et les variables
environnementales. On extrait les noms des colonnes de hull puis les
variables environnementales puis on analyse la correspondance avec
envfit. On affiche les résultats (pouvant inclure des stats pour chaque
variables environnementale). On fait ensuite un graphique pour les
visualiser.

``` r
data.frame(names(hull))
```

    ##           names.hull.
    ## 1              Axis.1
    ## 2              Axis.2
    ## 3     sample.SampName
    ## 4          sample.Geo
    ## 5  sample.Description
    ## 6       sample.groupe
    ## 7         sample.Pres
    ## 8      sample.PicoEuk
    ## 9        sample.Synec
    ## 10   sample.Prochloro
    ## 11     sample.NanoEuk
    ## 12      sample.Crypto
    ## 13       sample.SiOH4
    ## 14         sample.NO2
    ## 15         sample.NO3
    ## 16         sample.NH4
    ## 17         sample.PO4
    ## 18          sample.NT
    ## 19          sample.PT
    ## 20        sample.Chla
    ## 21           sample.T
    ## 22           sample.S
    ## 23     sample.Sigma_t

``` r
env <- hull[, 13:23]
ef <- vegan::envfit(physeq_clr_nmds, env, permu = 1000)
ef
```

    ## 
    ## ***VECTORS
    ## 
    ##                   NMDS1    NMDS2     r2   Pr(>r)    
    ## sample.SiOH4   -0.95409 -0.29952 0.2717 0.102897    
    ## sample.NO2     -0.44259 -0.89672 0.3271 0.062937 .  
    ## sample.NO3      0.94086  0.33880 0.2986 0.069930 .  
    ## sample.NH4     -0.48808 -0.87280 0.4484 0.018981 *  
    ## sample.PO4     -0.67398 -0.73875 0.2498 0.099900 .  
    ## sample.NT       0.02371 -0.99972 0.0526 0.688312    
    ## sample.PT      -0.61900 -0.78539 0.3745 0.035964 *  
    ## sample.Chla    -0.96843 -0.24930 0.2016 0.192807    
    ## sample.T       -0.87263 -0.48838 0.3250 0.051948 .  
    ## sample.S       -0.93218 -0.36199 0.7607 0.000999 ***
    ## sample.Sigma_t -0.96163 -0.27437 0.2116 0.205794    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 1000

``` r
plot(physeq_clr_nmds, type = "t", display = "sites")
plot(ef, p.max = 0.05)
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-53-1.png)<!-- -->

7.  On fait un test PERMANOVA pour évaluer si la variable Geo a un effet
    significatif sur la structure des communautés microbiennes
    représentée par les distances dans physeq_clr_dist. Il utilise
    ensuite les résultats de cette analyse pour déterminer si les
    différences observées entre les groupes définis par la variable Geo
    sont statistiquement significatives.

``` r
metadata <- data.frame(sample_data(physeq_clr))
results_permanova <- vegan::adonis2(physeq_clr_dist ~ Geo,
                                    data = metadata,
                                    perm = 1000)
results_permanova
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## vegan::adonis2(formula = physeq_clr_dist ~ Geo, data = metadata, permutations = 1000)
    ##          Df SumOfSqs      R2      F   Pr(>F)   
    ## Geo       1   1135.5 0.20329 4.0825 0.001998 **
    ## Residual 16   4450.1 0.79671                   
    ## Total    17   5585.6 1.00000                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

On évalue l’homogénéité des dispersions de groupes définis par la
variable Geo pour déterminer si la variabilité des distances entre les
échantillons au sein des groupes est similaire. Cela est fait avec la
fonction betadisper suivie d’une ANOVA. Ensuite, on réalise une analyse
PERMANOVA sur les données transformées physeq_clr_asv pour tester si la
structure des communautés microbiennes diffère significativement entre
les groupes définis par Geo. Après, on identifie et on affiche les 10
ASV les plus contributifs à ces différences, sous forme d’un graphique à
barres horizontales.

``` r
anova(vegan::betadisper(physeq_clr_dist, metadata$Geo))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## Groups     1 49.657  49.657  13.915 0.001822 **
    ## Residuals 16 57.096   3.569                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
permanova <- vegan::adonis(t(physeq_clr_asv) ~ Geo,
                            data = metadata,
                            permutations = 1000,
                            method = "euclidean")
```

    ## 'adonis' will be deprecated: use 'adonis2' instead

``` r
coef <- coefficients(permanova)["Geo1",]

top.coef <- coef[rev(order(abs(coef)))[1:10]]

par(mar = c(3, 14, 2, 1))

barplot(sort(top.coef),
        horiz = TRUE,
        las = 1,
        main = "Top taxa",
        cex.names = 0.7)
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-55-1.png)<!-- -->

On fait pareil avec le S.

``` r
permanova_S <- vegan::adonis2(physeq_clr_dist ~ S,
                              data = metadata,
                              perm = 1000)
permanova_S
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## vegan::adonis2(formula = physeq_clr_dist ~ S, data = metadata, permutations = 1000)
    ##          Df SumOfSqs      R2      F   Pr(>F)    
    ## S         1   1294.1 0.23168 4.8247 0.000999 ***
    ## Residual 16   4291.5 0.76832                    
    ## Total    17   5585.6 1.00000                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Puis pareil avec le NH4.

``` r
permanova_NH4 <- vegan::adonis2(physeq_clr_dist ~ NH4,
                                data = metadata,
                                perm = 1000)
permanova_NH4
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## vegan::adonis2(formula = physeq_clr_dist ~ NH4, data = metadata, permutations = 1000)
    ##          Df SumOfSqs      R2      F  Pr(>F)  
    ## NH4       1    769.8 0.13782 2.5575 0.01598 *
    ## Residual 16   4815.8 0.86218                 
    ## Total    17   5585.6 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Puis pareil avec le PT.

``` r
permanova_PT <- vegan::adonis2(physeq_clr_dist ~ PT,
                               data = metadata,
                               perm = 1000)
permanova_PT
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## vegan::adonis2(formula = physeq_clr_dist ~ PT, data = metadata, permutations = 1000)
    ##          Df SumOfSqs      R2      F  Pr(>F)  
    ## PT        1    697.3 0.12483 2.2822 0.01898 *
    ## Residual 16   4888.3 0.87517                 
    ## Total    17   5585.6 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

On regarde ensuite les résultats par rapport à toutes les variables.

``` r
permanova_all <- vegan::adonis2(physeq_clr_dist ~ SiOH4 + NO2 + NO3 + NH4 + PO4 + NT + PT + Chla + T + S + Sigma_t,
                                by="margin",
                                data=metadata,
                                perm=1000)

permanova_all
```

    ## Permutation test for adonis under reduced model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## vegan::adonis2(formula = physeq_clr_dist ~ SiOH4 + NO2 + NO3 + NH4 + PO4 + NT + PT + Chla + T + S + Sigma_t, data = metadata, permutations = 1000, by = "margin")
    ##          Df SumOfSqs      R2      F Pr(>F)
    ## SiOH4     1    291.5 0.05219 1.0594 0.3217
    ## NO2       1    243.4 0.04357 0.8846 0.5425
    ## NO3       1    239.8 0.04293 0.8715 0.5634
    ## NH4       1    253.2 0.04533 0.9203 0.5145
    ## PO4       1    232.7 0.04165 0.8456 0.5914
    ## NT        1    234.6 0.04200 0.8527 0.6034
    ## PT        1    234.9 0.04206 0.8539 0.5734
    ## Chla      1    200.8 0.03594 0.7296 0.7692
    ## T         1    285.9 0.05118 1.0390 0.3856
    ## S         1    286.2 0.05124 1.0402 0.3816
    ## Sigma_t   1    285.3 0.05108 1.0370 0.3896
    ## Residual  6   1650.8 0.29555              
    ## Total    17   5585.6 1.00000

On calcule la matrice de corrélation de Spearman pour certains colonnes
de metadata. On évalue ensuite la significativité de ces corrélations en
effectuant des tests stats. On visualise la matrice de corrélation en
utilisant un “corrplot”, où seules les corrélations significatives sont
affichées. Les corrélations non significatives sont laissées en blanc
dans le graphique.

``` r
cor_metadadata <- cor(metadata[, 11:21], method = "spearman")

cor_mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p_mat <- matrix(NA, n, n)
  diag(p_mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = "spearman", ...)
      p_mat[i, j] <- p_mat[j, i] <- tmp$p.value
    }
  }
  colnames(p_mat) <- rownames(p_mat) <- colnames(mat)
  p_mat
}

p_mat <- cor_mtest(metadata[, 11:21])
```

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(mat[, i], mat[, j], method = "spearman", ...):
    ## Cannot compute exact p-value with ties

``` r
corrplot::corrplot(cor_metadadata,
                   type = "upper",
                   order = "hclust",
                   p.mat = p_mat,
                   sig.level = 0.05,
                   insig = "blank")
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-60-1.png)<!-- -->

On effectue une analyse PERMANOVA pour évaluer l’effet de plusieurs
variables environnementales (S, NO3, NT, Chla, T) sur la structure des
communautés microbiennes représentée par les distances dans
physeq_clr_dist. L’argument by = “margin” signifie que l’effet de chaque
variable est évalué en tenant compte des autres. Les résultats de cette
analyse sont ensuite affichés pour déterminer si les différences
observées associées à ces variables environnementales sont
statistiquement significatives.

``` r
permanova_cor_pars <- vegan::adonis2(physeq_clr_dist ~ S + NO3 + NT + Chla + T,
                                     by = "margin",
                                     data = metadata,
                                     perm = 1000)
permanova_cor_pars
```

    ## Permutation test for adonis under reduced model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 1000
    ## 
    ## vegan::adonis2(formula = physeq_clr_dist ~ S + NO3 + NT + Chla + T, data = metadata, permutations = 1000, by = "margin")
    ##          Df SumOfSqs      R2      F  Pr(>F)  
    ## S         1    568.2 0.10173 2.0376 0.03896 *
    ## NO3       1    215.5 0.03858 0.7727 0.70230  
    ## NT        1    263.4 0.04716 0.9446 0.46653  
    ## Chla      1    170.9 0.03060 0.6129 0.94605  
    ## T         1    304.5 0.05452 1.0921 0.31169  
    ## Residual 12   3346.3 0.59910                 
    ## Total    17   5585.6 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

On teste si les différences entre les groupes définis par la variable
Geo sont plus grandes que ce que l’on pourrait attendre par hasard,
étant donné les distances ou dissimilarités entre les échantillons dans
physeq_clr_dist. Si le résultat est significatif, cela suggère que la
variable Geo a un effet sur la structure des communautés microbiennes.

``` r
vegan::anosim(physeq_clr_dist, metadata$Geo, permutations = 1000)
```

    ## 
    ## Call:
    ## vegan::anosim(x = physeq_clr_dist, grouping = metadata$Geo, permutations = 1000) 
    ## Dissimilarity: euclidean 
    ## 
    ## ANOSIM statistic R: 0.5628 
    ##       Significance: 0.001998 
    ## 
    ## Permutation: free
    ## Number of permutations: 1000

8.  On effectue une analyse de la redondance pour explorer les relations
    entre la structure des communautés microbiennes dans physeq_clr_asv
    et un ensemble de variables environnementales contenues dans
    certaines colonnes de metadata. On affiche ensuite un résumé des
    premiers résultats de cette analyse.

``` r
spe_rda <- vegan::rda(t(physeq_clr_asv) ~ .,
                      metadata[, 11:21])
head(summary(spe_rda))
```

    ## 
    ## Call:
    ## rda(formula = t(physeq_clr_asv) ~ SiOH4 + NO2 + NO3 + NH4 + PO4 +      NT + PT + Chla + T + S + Sigma_t, data = metadata[, 11:21]) 
    ## 
    ## Partitioning of variance:
    ##               Inertia Proportion
    ## Total          328.56     1.0000
    ## Constrained    231.46     0.7044
    ## Unconstrained   97.11     0.2956
    ## 
    ## Eigenvalues, and their contribution to the variance 
    ## 
    ## Importance of components:
    ##                          RDA1     RDA2     RDA3     RDA4     RDA5     RDA6
    ## Eigenvalue            85.2928 30.29173 20.29415 18.85659 15.83909 12.98651
    ## Proportion Explained   0.2596  0.09219  0.06177  0.05739  0.04821  0.03952
    ## Cumulative Proportion  0.2596  0.35179  0.41355  0.47094  0.51915  0.55868
    ##                           RDA7     RDA8     RDA9   RDA10   RDA11      PC1
    ## Eigenvalue            11.78027 10.97738 10.18119 7.94385 7.01222 28.88564
    ## Proportion Explained   0.03585  0.03341  0.03099 0.02418 0.02134  0.08791
    ## Cumulative Proportion  0.59453  0.62794  0.65893 0.68310 0.70445  0.79236
    ##                            PC2     PC3      PC4      PC5     PC6
    ## Eigenvalue            16.45693 16.3958 15.58129 11.19715 8.59184
    ## Proportion Explained   0.05009  0.0499  0.04742  0.03408 0.02615
    ## Cumulative Proportion  0.84245  0.8923  0.93977  0.97385 1.00000
    ## 
    ## Accumulated constrained eigenvalues
    ## Importance of components:
    ##                          RDA1    RDA2     RDA3     RDA4     RDA5     RDA6
    ## Eigenvalue            85.2928 30.2917 20.29415 18.85659 15.83909 12.98651
    ## Proportion Explained   0.3685  0.1309  0.08768  0.08147  0.06843  0.05611
    ## Cumulative Proportion  0.3685  0.4994  0.58706  0.66853  0.73696  0.79307
    ##                          RDA7     RDA8     RDA9   RDA10  RDA11
    ## Eigenvalue            11.7803 10.97738 10.18119 7.94385 7.0122
    ## Proportion Explained   0.0509  0.04743  0.04399 0.03432 0.0303
    ## Cumulative Proportion  0.8440  0.89139  0.93538 0.96970 1.0000
    ## 
    ## Scaling 2 for species and site scores
    ## * Species are scaled proportional to eigenvalues
    ## * Sites are unscaled: weighted dispersion equal on all dimensions
    ## * General scaling constant of scores:  8.645047 
    ## 
    ## 
    ## Species scores
    ## 
    ##                                                  RDA1      RDA2     RDA3
    ## ASV1_Cyanobiaceae_Synechococcus CC9902        -0.1033  0.108773  0.04666
    ## ASV2_Pseudoalteromonadaceae_Pseudoalteromonas -0.7807 -0.229145 -0.22860
    ## ASV3_Clade I_Clade Ia                         -0.2568  0.002182 -0.22536
    ## ASV4_NA_NA                                    -0.6996  0.193071  0.23547
    ## ASV5_Clade I_Clade Ia                          0.5264 -0.195773  0.23032
    ## ASV6_Clade II_NA                              -0.2542 -0.344583 -0.32380
    ## ....                                                                    
    ##                                                   RDA4     RDA5     RDA6
    ## ASV1_Cyanobiaceae_Synechococcus CC9902         0.12535 -0.01552  0.06487
    ## ASV2_Pseudoalteromonadaceae_Pseudoalteromonas -0.33352  0.13369  0.08880
    ## ASV3_Clade I_Clade Ia                          0.04191 -0.04528 -0.11436
    ## ASV4_NA_NA                                    -0.20648 -0.23531  0.06807
    ## ASV5_Clade I_Clade Ia                          0.05792  0.40196 -0.22286
    ## ASV6_Clade II_NA                               0.31352 -0.10920 -0.06137
    ## ....                                                                    
    ## 
    ## 
    ## Site scores (weighted sums of species scores)
    ## 
    ##        RDA1     RDA2    RDA3    RDA4     RDA5    RDA6
    ## S11B -1.703 -1.23820  2.9437  0.2362  1.13728 -0.4405
    ## S1B   2.565 -0.13340 -0.7868  5.7453  3.30268 -3.3657
    ## S2B   3.022 -2.96571  0.4021  0.9802 -3.09213  0.9282
    ## S2S  -1.731 -1.82618  2.0707  0.3281 -0.66853 -1.6638
    ## S3B   3.624 -1.55655 -1.2829  2.0701 -2.02586  1.7347
    ## S3S   3.165 -0.08923  2.8998 -2.0441 -0.08464  2.0314
    ## ....                                                 
    ## 
    ## 
    ## Site constraints (linear combinations of constraining variables)
    ## 
    ##         RDA1    RDA2    RDA3    RDA4    RDA5    RDA6
    ## S11B -1.2105 -0.7764  3.0649  0.2199  1.2569  0.7586
    ## S1B   1.7387  0.3983 -0.3817  5.4943  3.2411 -2.7484
    ## S2B   2.0536 -3.3237  0.6260  1.4897 -2.8936  0.1774
    ## S2S   0.5936 -2.0609  1.1588  0.1736 -0.8183 -1.8069
    ## S3B   4.1498 -1.1569 -1.6837  1.1942 -2.4216  2.5295
    ## S3S   2.0704 -0.1285  3.6947 -1.1733  0.3885  1.8438
    ## ....                                                
    ## 
    ## 
    ## Biplot scores for constraining variables
    ## 
    ##             RDA1     RDA2     RDA3     RDA4     RDA5      RDA6
    ## SiOH4   -0.57424 -0.21106 -0.25450 -0.25678 -0.02349 -0.213981
    ## NO2     -0.51463 -0.10086 -0.08171  0.34294  0.35340  0.013696
    ## NO3      0.59878  0.05632 -0.04267 -0.02065 -0.30772  0.095439
    ## NH4     -0.63097 -0.49073 -0.01146 -0.07457  0.25646  0.259440
    ## PO4     -0.49369 -0.05367 -0.31521  0.04459  0.19877  0.304690
    ## NT       0.02778 -0.05873 -0.28198  0.59590  0.14825 -0.392684
    ## PT      -0.61634 -0.27995 -0.01129  0.12013  0.07328 -0.533916
    ## Chla    -0.47936 -0.07832 -0.06090 -0.01293 -0.11376  0.179421
    ## T       -0.57485  0.21879  0.26190  0.53662 -0.42902  0.007286
    ## S       -0.93622  0.00815 -0.06712  0.05543  0.04078  0.183950
    ## Sigma_t -0.52380 -0.20293 -0.31121 -0.40702  0.43162  0.205711

On calcule et on affiche le coefficient de détermination R2 et R2 ajusté
de l’analyse réalisée juste avant. Ces coefficients indiquent la
proportion de la variance dans les données de la communauté
microbiologique qui est expliquée par les variables environnementales
utilisées dans la RDA.

``` r
R2 <- vegan::RsquareAdj(spe_rda)$r.squared
R2
```

    ## [1] 0.7044457

``` r
R2adj <- vegan::RsquareAdj(spe_rda)$adj.r.squared
R2adj
```

    ## [1] 0.1625961

On réalise des tests de permutation ANOVA sur les résultats de l’analyse
de la redondance pour évaluer la significativité des axes et de
l’ensemble du modèle. On effectue d’abord un test sur l’ensemble du
modèle RDA, puis un test pour chaque axe de la RDA individuellement. Ces
tests utilisent des permutations pour évaluer la significativité.

``` r
anova(spe_rda, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = t(physeq_clr_asv) ~ SiOH4 + NO2 + NO3 + NH4 + PO4 + NT + PT + Chla + T + S + Sigma_t, data = metadata[, 11:21])
    ##          Df Variance      F Pr(>F)  
    ## Model    11  231.456 1.3001  0.096 .
    ## Residual  6   97.109                
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(spe_rda, by = "axis", step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = t(physeq_clr_asv) ~ SiOH4 + NO2 + NO3 + NH4 + PO4 + NT + PT + Chla + T + S + Sigma_t, data = metadata[, 11:21])
    ##          Df Variance      F Pr(>F)
    ## RDA1      1   85.293 5.2699  0.108
    ## RDA2      1   30.292 1.8716  0.614
    ## RDA3      1   20.294 1.2539  0.999
    ## RDA4      1   18.857 1.1651  1.000
    ## RDA5      1   15.839 0.9786  1.000
    ## RDA6      1   12.987 0.8024  1.000
    ## RDA7      1   11.780 0.7279  1.000
    ## RDA8      1   10.977 0.6783  1.000
    ## RDA9      1   10.181 0.6291  0.999
    ## RDA10     1    7.944 0.4908  0.993
    ## RDA11     1    7.012 0.4333  0.911
    ## Residual  6   97.109

On calcule les facteurs d’inflation de la variance (VIF) pour évaluer la
multicollinéarité entre les variables environnementales utilisées dans
la RDA et on réalise une sélection séquentielle de variables (stepwise)
en utilisant une méthode orientée vers l’avant (forward) pour déterminer
quelles variables environnementales contribuent de manière significative
à la structure des communautés microbiennes.

``` r
vegan::vif.cca(spe_rda)
```

    ##       SiOH4         NO2         NO3         NH4         PO4          NT 
    ##    4.066588    3.489186    3.634643   16.867288    8.819736    4.908553 
    ##          PT        Chla           T           S     Sigma_t 
    ##    6.835572    2.264012 5417.455601 8388.550079 6878.896122

``` r
step_forward <- vegan::ordiR2step(vegan::rda(t(physeq_clr_asv) ~ 1,
                                             data = metadata[, 11:21]),
                                  scope = formula(spe_rda),
                                  direction = "forward",
                                  pstep = 1000)
```

    ## Step: R2.adj= 0 
    ## Call: t(physeq_clr_asv) ~ 1 
    ##  
    ##                 R2.adjusted
    ## + S              0.18366030
    ## <All variables>  0.16259613
    ## + NH4            0.08392874
    ## + PT             0.07013415
    ## + T              0.06719602
    ## + NO3            0.05904665
    ## + SiOH4          0.05787026
    ## + Sigma_t        0.05002017
    ## + NO2            0.03846019
    ## + PO4            0.03190148
    ## + Chla           0.02451726
    ## <none>           0.00000000
    ## + NT            -0.01448677

On se concentre ensuite plus particulièrement sur la variable S pour
explorer sa relation avec la structure des communautés microbiennes dans
physeq_clr_asv. Ensuite, on effectue des tests de permutation ANOVA sur
les résultats de cette RDA pour évaluer la significativité de l’ensemble
du modèle ainsi que de chaque axe individuellement.

``` r
spe_rda_pars <- vegan::rda(t(physeq_clr_asv) ~ S, data = metadata[, 11:21])
anova(spe_rda_pars, step = 1000)
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = t(physeq_clr_asv) ~ S, data = metadata[, 11:21])
    ##          Df Variance      F Pr(>F)    
    ## Model     1   76.122 4.8247  0.001 ***
    ## Residual 16  252.443                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(spe_rda_pars, step = 1000, by = "axis")
```

    ## Permutation test for rda under reduced model
    ## Forward tests for axes
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: rda(formula = t(physeq_clr_asv) ~ S, data = metadata[, 11:21])
    ##          Df Variance      F Pr(>F)    
    ## RDA1      1   76.122 4.8247  0.001 ***
    ## Residual 16  252.443                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

On calcule stocke le coefficient de détermination R2 et R2 ajusté pour
les résultats de l’Analyse de la Redondance (RDA) basée sur la variable
“S”. Puis on calcule les facteurs d’inflation de la variance (VIF) pour
évaluer la multicollinéarité entre les variables environnementales de
l’analyse RDA originale. Enfin, on calcule les VIF pour l’analyse RDA
basée uniquement sur la variable “S” pour évaluer la multicollinéarité.

``` r
R2adj_pars <- vegan::RsquareAdj(spe_rda_pars)$adj.r.squared
vegan::vif.cca(spe_rda)
```

    ##       SiOH4         NO2         NO3         NH4         PO4          NT 
    ##    4.066588    3.489186    3.634643   16.867288    8.819736    4.908553 
    ##          PT        Chla           T           S     Sigma_t 
    ##    6.835572    2.264012 5417.455601 8388.550079 6878.896122

``` r
vegan::vif.cca(spe_rda_pars)
```

    ## S 
    ## 1

On fait un plot qui permet de visualiser la position des échantillons
(sites) dans l’espace RDA, distingués par la variable “Geo”, les six ASV
les plus contributives à la variation observée et l’influence de la
salinité sur la structure des communautés microbiennes. Le graphique
résultant aide à interpréter comment les échantillons se regroupent en
fonction de la variable Geo et à comprendre l’effet des espèces les plus
contributives et de la salinité sur cette structure.

``` r
ii <- summary(spe_rda_pars)

sp <- as.data.frame(ii$species[, 1:2]) * 2
sp_top <- sp[order(abs(sp$RDA1), decreasing = TRUE), ][1:6, ]

st <- as.data.frame(ii$sites[, 1:2])
st <- merge(st,
      metadata["Geo"],
      by = "row.names")

yz <- t(as.data.frame(ii$biplot[, 1:2]))
row.names(yz) <- "Salinity"
yz <- as.data.frame(yz)

eigen_values <- format(100 *ii$cont[[1]][2,], digits=4)

ggplot() +
  geom_point(data = st, size = 4,
             aes(x = RDA1, y = PC1,
                 shape = Geo, fill = Geo)) +
  scale_shape_manual(values = c(21:25)) +
  geom_segment(data = sp_top,
               arrow = arrow(angle = 22.5,
                             length = unit(0.35, "cm"),
                             type = "closed"),
               linetype = 1, size = 0.6, colour = "red",
               aes(x = 0, y = 0, xend = RDA1, yend = PC1)) +
  ggrepel::geom_text_repel(data = sp_top,
                           aes(x = RDA1, y = PC1, label = row.names(sp_top))) +
  geom_segment(data = yz,
               arrow = arrow(angle = 22.5,
                             length = unit(0.35,"cm"),
                             type = "closed"),
               linetype = 1, size = 0.6, colour = "blue",
               aes(x = 0, y = 0, xend = RDA1, yend = PC1)) +
  ggrepel::geom_text_repel(data = yz, aes(RDA1, PC1, label=row.names(yz)))+
  labs(x = paste("RDA 1 (", eigen_values[1], "%)", sep = ""),
       y = paste("PC 1 (", eigen_values[2], "%)", sep = ""))+
  geom_hline(yintercept = 0,linetype = 3,size = 1) + 
  geom_vline(xintercept = 0,linetype = 3,size = 1)+
  guides(shape = guide_legend(title = NULL,
         color = "black"),
         fill = guide_legend(title = NULL))+
  theme_bw() +
  theme(panel.grid = element_blank())
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-69-1.png)<!-- -->

On indique à R de copier le dossier beta_diversity dans le /data du menu
home.

``` bash
cp -R /home/rstudio/ADM2023_tutoriel/course-material-main/data/beta_diversity ./data
```

On lit le fichier spatial_distance.rds et on calcule une matrice de
distance à partir de cet objet.

``` r
ANF_km <- readRDS(here::here("data","beta_diversity","spatial_distance.rds"))
ANF_km_dist <- dist(ANF_km)
```

On effectue une analyse de décroissance de la dissimilarité. L’idée est
d’examiner comment la dissimilarité écologique (ici basée sur la
transformation CLR) entre les échantillons change en fonction de leur
distance spatiale. On utilise le modèle exponentiel de décroissance pour
ajuster les données. On affiche un graphique de dispersion des distances
spatiales par rapport à la dissimilarité. On superpose le modèle de
décroissance exponentiel ajusté sur le graphique de dispersion. On met
une légende indiquant les paramètres et statistiques du modèle ajusté.
L’objectif est de voir si, à mesure que les échantillons sont plus
éloignés les uns des autres dans l’espace, ils deviennent également plus
dissimilaires écologiquement.

``` r
ANF_decay_exp <- betapart::decay.model(physeq_clr_dist/100,
                                       ANF_km_dist,
                                       y.type="dissim",
                                       model.type="exp",
                                       perm=100)
plot(ANF_km_dist, physeq_clr_dist/100,
     ylim=c(0, max(physeq_clr_dist/100)),
     xlim=c(0, max(ANF_km_dist)),
     xlab = "Distance (km)", ylab = "Dissimilarity (CLR)")

betapart::plot.decay(ANF_decay_exp, col = "blue",
                     remove.dots = TRUE, add = TRUE)

legend("bottomright",
       paste("exp: (Beta =", round(ANF_decay_exp$second.parameter, 4),
             ", Rsqr =", round(ANF_decay_exp$pseudo.r.squared, 2),
             ", p =", round(ANF_decay_exp$p.value, 2)),
       fill = "blue")
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-72-1.png)<!-- -->

On fait une analyse de Mantel régressée par morceaux. On calcule une
matrice de distance euclidienne pour les données physeq_clr. On calcule
également une matrice de distance pour les distances spatiales ANF_km de
la même manière. Puis, on calcule une matrice de distance pour certaines
variables environnementales présentes dans metadata. On utilise ensuite
ces matrices pour effectuer une analyse de Mantel afin d’évaluer comment
la structure des communautés microbiennes (distance physeq_clr) est
expliquée par les distances spatiales et les variables environnementales
simultanément. L’objectif est de déterminer dans quelle mesure les
variables environnementales et spatiales influencent conjointement la
structure des communautés microbiennes.

``` r
physeq_clr_dist_square <- phyloseq::distance(physeq_clr,
                                             method = "euclidean",
                                             diag = TRUE,
                                             upper = TRUE)
ANF_km_dist_square <- dist(ANF_km, diag = TRUE, upper = TRUE)
envdata <- dist(metadata[,11:21], diag = TRUE, upper = TRUE)
ecodist::MRM(physeq_clr_dist_square ~ envdata + ANF_km_dist_square, nperm=1000)
```

    ## $coef
    ##                    physeq_clr_dist_square  pval
    ## Int                           19.45167946 0.908
    ## envdata                        1.28567618 0.001
    ## ANF_km_dist_square             0.01828172 0.001
    ## 
    ## $r.squared
    ##       R2     pval 
    ## 0.366774 0.001000 
    ## 
    ## $F.test
    ##        F   F.pval 
    ## 43.44112  0.00100

On effectue une analyse de Mantel régressée par morceaux (MRM) pour
évaluer la relation entre la structure des communautés microbiennes
(distance physeq_clr) et les variables environnementales. Puis on
réalise une autre analyse MRM pour évaluer la relation entre la
structure des communautés microbiennes (distance physeq_clr) et les
distances spatiales (ANF_km). Enfin, on utilise la fonction varPart du
package modEvA pour effectuer une partition de la variance, afin de
décomposer la contribution relative des effets environnementaux et de la
limitation de la dispersion (probablement en termes d’influence
spatiale) sur la structure des communautés microbiennes. L’objectif est
de décomposer l’influence relative des facteurs environnementaux et
spatiaux sur la structure des communautés microbiennes.

``` r
ecodist::MRM(physeq_clr_dist_square ~ envdata, nperm=1000)
```

    ## $coef
    ##         physeq_clr_dist_square  pval
    ## Int                  21.042622 0.967
    ## envdata               1.609333 0.001
    ## 
    ## $r.squared
    ##        R2      pval 
    ## 0.2122659 0.0010000 
    ## 
    ## $F.test
    ##        F   F.pval 
    ## 40.68905  0.00100

``` r
ecodist::MRM(physeq_clr_dist_square ~ ANF_km_dist_square, nperm=1000)
```

    ## $coef
    ##                    physeq_clr_dist_square  pval
    ## Int                           22.34249373 0.354
    ## ANF_km_dist_square             0.02210456 0.004
    ## 
    ## $r.squared
    ##        R2      pval 
    ## 0.2384328 0.0040000 
    ## 
    ## $F.test
    ##        F   F.pval 
    ## 47.27535  0.00400

``` r
modEvA::varPart(A = 0.212, B = 0.238, AB = 0.366,
                A.name = "Environmental",
                B.name = "Dispersal limitation")
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-74-1.png)<!-- -->

    ##                                    Proportion
    ## Environmental                           0.128
    ## Dispersal limitation                    0.154
    ## Environmental_Dispersal limitation      0.084
    ## Unexplained                             0.634

9.  On réalise une analyse LEfSe sur les données phylogénétiques
    fournies par physeq. La fonction run_lefse du package
    microbiomeMarker est utilisée pour identifier les ASV qui sont
    statistiquement différentes entre les groupes définis par la
    variable “Geo”. La normalisation est faite en utilisant la méthode
    “Counts Per Million” (CPM). Des seuils spécifiques sont définis pour
    les tests de Wilcoxon et Kruskal-Wallis. L’analyse est configurée
    pour identifier des marqueurs dans des groupes multiples et un seuil
    pour la taille de l’effet est défini. Les résultats de l’analyse
    LEfSe sont ensuite convertis en une dataframe pour faciliter
    l’affichage et la manipulation. L’objectif est d’identifier quels
    microorganismes (ou groupes de microorganismes) sont
    significativement associés à chacun des groupes définis par la
    variable “Geo”.

``` r
mm_lefse <- microbiomeMarker::run_lefse(physeq, norm = "CPM",
                                        wilcoxon_cutoff = 0.01,
                                        group = "Geo",
                                        taxa_rank = "none",
                                        kw_cutoff = 0.01,
                                        multigrp_strat = TRUE,
                                        lda_cutoff = 4)
```

    ## Registered S3 method overwritten by 'httr':
    ##   method         from  
    ##   print.response rmutil

    ## Registered S3 methods overwritten by 'treeio':
    ##   method              from    
    ##   MRCA.phylo          tidytree
    ##   MRCA.treedata       tidytree
    ##   Nnode.treedata      tidytree
    ##   Ntip.treedata       tidytree
    ##   ancestor.phylo      tidytree
    ##   ancestor.treedata   tidytree
    ##   child.phylo         tidytree
    ##   child.treedata      tidytree
    ##   full_join.phylo     tidytree
    ##   full_join.treedata  tidytree
    ##   groupClade.phylo    tidytree
    ##   groupClade.treedata tidytree
    ##   groupOTU.phylo      tidytree
    ##   groupOTU.treedata   tidytree
    ##   inner_join.phylo    tidytree
    ##   inner_join.treedata tidytree
    ##   is.rooted.treedata  tidytree
    ##   nodeid.phylo        tidytree
    ##   nodeid.treedata     tidytree
    ##   nodelab.phylo       tidytree
    ##   nodelab.treedata    tidytree
    ##   offspring.phylo     tidytree
    ##   offspring.treedata  tidytree
    ##   parent.phylo        tidytree
    ##   parent.treedata     tidytree
    ##   root.treedata       tidytree
    ##   rootnode.phylo      tidytree
    ##   sibling.phylo       tidytree

    ## Registered S3 method overwritten by 'gplots':
    ##   method         from     
    ##   reorder.factor DescTools

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

``` r
mm_lefse_table <- data.frame(mm_lefse@marker_table)
mm_lefse_table
```

    ##          feature enrich_group   ef_lda       pvalue         padj
    ## marker1    ASV11        North 4.746047 0.0015574784 0.0015574784
    ## marker2    ASV12        North 4.699727 0.0045142882 0.0045142882
    ## marker3    ASV10        North 4.661376 0.0022950748 0.0022950748
    ## marker4    ASV18        North 4.460565 0.0045142882 0.0045142882
    ## marker5    ASV35        North 4.183570 0.0045142882 0.0045142882
    ## marker6    ASV49        North 4.025863 0.0045142882 0.0045142882
    ## marker7     ASV2        South 4.950924 0.0039173223 0.0039173223
    ## marker8     ASV8        South 4.706448 0.0020814438 0.0020814438
    ## marker9     ASV7        South 4.670957 0.0010275895 0.0010275895
    ## marker10    ASV3        South 4.433849 0.0091897421 0.0091897421
    ## marker11   ASV13        South 4.406032 0.0073724319 0.0073724319
    ## marker12   ASV27        South 4.333577 0.0008112059 0.0008112059

On génère et on affiche deux types de graphiques basés sur les résultats
de l’analyse LEfSe réalisée précédemment sur les données
phylogénétiques. plot_ef_bar crée un graphique à barres montrant la
taille de l’effet (LDA score) pour chaque ASV identifié comme étant
significativement différent entre les groupes. Plus la barre est haute
(ou basse), plus la différence est marquée. plot_abundance génère un
graphique à barres empilées montrant l’abondance relative des
microorganismes identifiés comme marqueurs entre les différents groupes
définis par la variable “Geo”. Enfin, grid.arrange combine et affiche
ces deux graphiques côte à côte pour une comparaison visuelle.
L’objectif est de visualiser à la fois la taille de l’effet et
l’abondance des microorganismes identifiés comme étant significativement
différents entre les groupes.

``` r
p_LDAsc <- microbiomeMarker::plot_ef_bar(mm_lefse)
y_labs <- ggplot_build(p_LDAsc)$layout$panel_params[[1]]$y$get_labels()
p_abd <- microbiomeMarker::plot_abundance(mm_lefse, group = "Geo") +
  scale_y_discrete(limits = y_labs)
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
gridExtra::grid.arrange(p_LDAsc, p_abd, nrow = 1)
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-76-1.png)<!-- -->

On effectue une analyse ANCOM-BC sur les données phylogénétiques
fournies par physeq. La fonction run_ancombc_patched est utilisée pour
identifier les espèces (ou ASV/OTU) qui sont statistiquement différentes
entre les groupes définis par la variable “Geo”. Un seuil de p-valeur
est défini pour déterminer la significativité statistique. La correction
pour les tests multiples est réalisée en utilisant la méthode FDR (False
Discovery Rate). Les résultats de l’analyse ANCOM-BC sont ensuite
convertis en une dataframe pour faciliter l’affichage et la
manipulation. L’objectif est d’identifier quels microorganismes (ou
groupes de microorganismes) montrent des abondances significativement
différentes entre les groupes définis par la variable “Geo”, tout en
corrigeant pour les biais potentiels dans l’analyse des données
composées comme celles des données microbiennes.

``` r
mm_ancombc <- run_ancombc_patched(
  physeq,
  group = "Geo",
  taxa_rank = "none",
  pvalue_cutoff = 0.001,
  p_adjust = "fdr"
)
```

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## 'ancombc' has been fully evolved to 'ancombc2'. 
    ## Explore the enhanced capabilities of our refined method!

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## `tax_level` is not specified 
    ## No agglomeration will be performed
    ## Otherwise, please specify `tax_level` by one of the following:

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Warning: The group variable has < 3 categories 
    ## The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated

    ## Loading required package: foreach

    ## Loading required package: rngtools

``` r
mm_ancombc_table <- data.frame(mm_ancombc@marker_table)
mm_ancombc_table
```

    ##          feature enrich_group      ef_W       pvalue         padj
    ## marker1     ASV2        South  3.980197 6.885820e-05 7.230111e-04
    ## marker2     ASV7        South  4.341347 1.416118e-05 1.652137e-04
    ## marker3     ASV8        South  4.532481 5.829496e-06 1.020162e-04
    ## marker4    ASV10        North -4.775089 1.796277e-06 6.286968e-05
    ## marker5    ASV11        North -5.811580 6.188594e-09 3.249012e-07
    ## marker6    ASV12        North -4.466839 7.938375e-06 1.041912e-04
    ## marker7    ASV18        North -4.561024 5.090471e-06 1.020162e-04
    ## marker8    ASV27        South  5.874154 4.250091e-09 3.249012e-07
    ## marker9    ASV35        North -4.483869 7.330158e-06 1.041912e-04
    ## marker10   ASV49        North -4.680720 2.858686e-06 7.504051e-05

On affiche deux types de graphiques basés sur les résultats de l’analyse
ANCOM-BC réalisée précédemment sur les données phylogénétiques.
plot_ef_bar crée un graphique à barres qui illustre la taille de l’effet
pour chaque microorganisme (ou ASV/OTU) identifié comme étant
significativement différent entre les groupes. Les barres indiquent
l’importance de chaque différence détectée. plot_abundance produit un
graphique à barres empilées montrant l’abondance relative des
microorganismes identifiés comme marqueurs entre les différents groupes
définis par la variable “Geo”. Ensuite, grid.arrange combine et affiche
ces deux graphiques côte à côte pour une comparaison visuelle.
L’objectif est de visualiser à la fois la taille de l’effet et
l’abondance des microorganismes qui ont été identifiés comme étant
significativement différents entre les groupes à partir de l’analyse
ANCOM-BC.

``` r
an_ef <- microbiomeMarker::plot_ef_bar(mm_ancombc)
y_labs <- ggplot_build(an_ef)$layout$panel_params[[1]]$y$get_labels()
an_abd <- microbiomeMarker::plot_abundance(mm_ancombc, group = "Geo") +
  scale_y_discrete(limits = y_labs)
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
gridExtra::grid.arrange(an_ef, an_abd, nrow = 1)
```

![](ADM2023_tutoriel_beta_diversite_files/figure-gfm/unnamed-chunk-78-1.png)<!-- -->

On effectue une analyse ALDEx2 (ANalysis of Differential Abundance
taking into account the Library sIzE) sur les données phylogénétiques
fournies par physeq. La fonction run_aldex du package microbiomeMarker
est utilisée pour identifier les ASV qui sont statistiquement
différentes entre les groupes définis par la variable “Geo”. La
normalisation est effectuée en utilisant la méthode “Counts Per Million”
(CPM). La correction pour les tests multiples est réalisée en utilisant
la méthode FDR (False Discovery Rate). Les résultats de l’analyse ALDEx2
sont ensuite convertis en une dataframe pour faciliter l’affichage et la
manipulation. L’objectif est d’identifier quels microorganismes (ou
groupes de microorganismes) montrent des abondances significativement
différentes entre les groupes définis par la variable “Geo” en utilisant
l’approche ALDEx2, qui est spécialement conçue pour traiter les données
de séquençage à haut débit.

``` r
mm_aldex <- microbiomeMarker::run_aldex(physeq, group = "Geo",
                                        norm = "CPM",
                                        taxa_rank = "none",
                                        p_adjust = "fdr")
```

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'

    ## Also defined by 'tidytree'

    ## operating in serial mode

    ## Warning: Not all reads are integers, the reads are ceiled to integers.
    ##    Raw reads is recommended from the ALDEx2 paper.

    ## operating in serial mode

    ## computing center with all features

    ## New names:
    ## • `` -> `...1`
    ## • `` -> `...2`
    ## • `` -> `...3`
    ## • `` -> `...4`
    ## • `` -> `...5`
    ## • `` -> `...6`
    ## • `` -> `...7`
    ## • `` -> `...8`
    ## • `` -> `...9`
    ## • `` -> `...10`
    ## • `` -> `...11`
    ## • `` -> `...12`
    ## • `` -> `...13`
    ## • `` -> `...14`
    ## • `` -> `...15`
    ## • `` -> `...16`
    ## • `` -> `...17`
    ## • `` -> `...18`
    ## • `` -> `...19`
    ## • `` -> `...20`
    ## • `` -> `...21`
    ## • `` -> `...22`
    ## • `` -> `...23`
    ## • `` -> `...24`
    ## • `` -> `...25`
    ## • `` -> `...26`
    ## • `` -> `...27`
    ## • `` -> `...28`
    ## • `` -> `...29`
    ## • `` -> `...30`
    ## • `` -> `...31`
    ## • `` -> `...32`
    ## • `` -> `...33`
    ## • `` -> `...34`
    ## • `` -> `...35`
    ## • `` -> `...36`
    ## • `` -> `...37`
    ## • `` -> `...38`
    ## • `` -> `...39`
    ## • `` -> `...40`
    ## • `` -> `...41`
    ## • `` -> `...42`
    ## • `` -> `...43`
    ## • `` -> `...44`
    ## • `` -> `...45`
    ## • `` -> `...46`
    ## • `` -> `...47`
    ## • `` -> `...48`
    ## • `` -> `...49`
    ## • `` -> `...50`
    ## • `` -> `...51`
    ## • `` -> `...52`
    ## • `` -> `...53`
    ## • `` -> `...54`
    ## • `` -> `...55`
    ## • `` -> `...56`
    ## • `` -> `...57`
    ## • `` -> `...58`
    ## • `` -> `...59`
    ## • `` -> `...60`
    ## • `` -> `...61`
    ## • `` -> `...62`
    ## • `` -> `...63`
    ## • `` -> `...64`
    ## • `` -> `...65`
    ## • `` -> `...66`
    ## • `` -> `...67`
    ## • `` -> `...68`
    ## • `` -> `...69`
    ## • `` -> `...70`
    ## • `` -> `...71`
    ## • `` -> `...72`
    ## • `` -> `...73`
    ## • `` -> `...74`
    ## • `` -> `...75`
    ## • `` -> `...76`
    ## • `` -> `...77`
    ## • `` -> `...78`
    ## • `` -> `...79`
    ## • `` -> `...80`
    ## • `` -> `...81`
    ## • `` -> `...82`
    ## • `` -> `...83`
    ## • `` -> `...84`
    ## • `` -> `...85`
    ## • `` -> `...86`
    ## • `` -> `...87`
    ## • `` -> `...88`
    ## • `` -> `...89`
    ## • `` -> `...90`
    ## • `` -> `...91`
    ## • `` -> `...92`
    ## • `` -> `...93`
    ## • `` -> `...94`
    ## • `` -> `...95`
    ## • `` -> `...96`
    ## • `` -> `...97`
    ## • `` -> `...98`
    ## • `` -> `...99`
    ## • `` -> `...100`
    ## • `` -> `...101`
    ## • `` -> `...102`
    ## • `` -> `...103`
    ## • `` -> `...104`
    ## • `` -> `...105`
    ## • `` -> `...106`
    ## • `` -> `...107`
    ## • `` -> `...108`
    ## • `` -> `...109`
    ## • `` -> `...110`
    ## • `` -> `...111`
    ## • `` -> `...112`
    ## • `` -> `...113`
    ## • `` -> `...114`
    ## • `` -> `...115`
    ## • `` -> `...116`
    ## • `` -> `...117`
    ## • `` -> `...118`
    ## • `` -> `...119`
    ## • `` -> `...120`
    ## • `` -> `...121`
    ## • `` -> `...122`
    ## • `` -> `...123`
    ## • `` -> `...124`
    ## • `` -> `...125`
    ## • `` -> `...126`
    ## • `` -> `...127`
    ## • `` -> `...128`

``` r
mm_aldex_table <- data.frame(mm_aldex@marker_table)
mm_aldex_table
```

    ##         feature enrich_group ef_aldex       pvalue       padj
    ## marker1   ASV27        North 2.095543 0.0003582814 0.03277281
