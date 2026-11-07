# References — *Heliconius erato etylus*, *H. erato demophoon*, *H. himera*

Companion to `plot-erato-himera-map.R` / `erato-himera-distribution.pdf` and `localities-erato-himera.tsv`.
BibTeX: `references-erato-himera.bib`. All downloads/queries: 2026-10-21.

## 1. Occurrence data — sources, priority, curation

### Source priority (highest first)

1. **GBIF** — primary source. GBIF itself aggregates iNaturalist research-grade observations
   (dataset `50c9509d-22c7-4a22-a47d-8c48425ef4a7`): 86 of the 87 GBIF *himera* rows and 99 of
   100 sampled GBIF *demophoon* rows come from that dataset, so GBIF ⊇ most of iNat.
2. **iNaturalist** — supplement, **research-grade observations only** (`quality_grade = research`;
   `needs_id` excluded). An iNat record is added only when no TSV record lies within 0.005°,
   i.e. only records absent from GBIF enter the TSV (9 *demophoon*, 16 *himera*, 17 *etylus*).

### Counts and curation per taxon

| Taxon | GBIF (plotted) | iNaturalist research-grade (plotted) | Literature (TSV/docs only, not plotted) | Curation applied |
|---|---|---|---|---|
| *H. himera* | 87 | 16 | — | GBIF taxonKey 1900219, `hasCoordinate=true`, restricted to the dry-valley belt (−8.0 ≤ lat ≤ −3.5, −80.5 ≤ lon ≤ −77.0); out-of-range community records excluded (Colombia 3.52,−76.62; Pichincha −0.05/−0.18; Cuenca basin −3.26,−79.96). iNat taxon 790498, dedup vs GBIF. |
| *H. erato demophoon* | 111 | 9 | — | GBIF taxonKey 11869372, `hasCoordinate=true`, restricted to Costa Rica + Panama (misidentified records from Brazil, Mexico, Canada, Colombia excluded). iNat taxon 1552759, dedup vs GBIF. |
| *H. erato etylus* | **0** | 17 | 6 | **GBIF holds no *etylus* occurrences at all**: backbone subspecies keys 104430119 and 177943358 return 0 records (with or without coordinates), and GBIF's iNaturalist research-grade dataset contains 0 *etylus* records (vs 86 for *himera*) — iNat *etylus* identifications are not indexed by GBIF at subspecies level. iNat taxon 1136782: two verifiable outliers excluded as misidentifications (urban Quito, Parque La Carolina ~2,800 m, obs 63608868; Shushufindi Amazonian lowland = *lativitta* territory, obs 109304571). Literature: see §3. |

### Verification links

- GBIF, per record: `https://www.gbif.org/occurrence/<key>` (key in the TSV `source` column, `GBIF:<key>`). Reproducible downloads:
  - *H. himera*: https://www.gbif.org/occurrence/search?taxon_key=1900219&has_coordinate=true
  - *H. erato demophoon*: https://www.gbif.org/occurrence/search?taxon_key=11869372&has_coordinate=true&country=Costa+Rica&country=Panama
- iNaturalist, per record: `https://www.inaturalist.org/observations/<id>` (id in the TSV `source` column). Taxon maps: [etylus 1136782](https://www.inaturalist.org/taxa/1136782-Heliconius-erato-etylus), [demophoon 1552759](https://www.inaturalist.org/taxa/1552759-Heliconius-erato-demophoon), [himera 790498](https://www.inaturalist.org/taxa/790498-Heliconius-himera).
- Literature rows carry their citation in the TSV `source` column (see §5).

## 2. Genomes and sequenced-individual provenance

NCBI BioSample `geo_loc_name` / isolation source, accessed 2026-10-21 via the NCBI Datasets API.

| Taxon | Accession (assembly) | Isolate | geo_loc_name | Source |
|---|---|---|---|---|
| *H. erato etylus* | GCA_056362835.1 (Heety_Duke_1.0) | WOM-AMS-2026c | Panama: Gamboa | STRI captive colony; founders exported from Ecuador, permit PA-03-018-2025 (Mackay-Smith et al., 2026) |
| *H. erato demophoon* | GCA_056362945.2 (Hedem_Duke_1.1) | WOM-AMS-2026b | Panama: Gamboa | STRI captive colony, Panama-native, wild-infused, permit ARB-110-2022 (Mackay-Smith et al., 2026) |
| *H. himera* | GCA_013421025.1 (ASM1342102v1) | 6363 | Ecuador | Inbred lab line; Loja-province stock (Jiggins et al., 1996; Van Belleghem et al., 2021) |
| *H. himera* | GCA_013421065.1 (ASM1342106v1) | 6900 | Ecuador | Inbred lab line; as above |
| *H. erato* × *H. himera* F1 | GCA_013421005.1 (ASM1342100v1) | 6700 | Panama (lab) | Laboratory F1 |

## 3. Literature / collection localities (*etylus* only; not plotted)

1. Type locality Guadalquiza, Pastaza, Ecuador (Salvin, 1871; type locality per Funet index).
2. Llandia, route Tena–Puyo km 25, Napo (LepiWiki photo record, 2018).
3. Cab. Yankuam, Río Nangaritza, Zamora Chinchipe, 1,000–1,100 m (ETH Zürich coll. AZ23-0322/0324, 2023).
4. San Carlos de las Minas, Zamora Chinchipe, 1,400 m (ETH Zürich coll. AZ23-0323, 2023).
5. Río Zamora, Loja, 3°55′S 78°50′W — *etylus* stock of the *erato*×*himera* mapping cross (Kapan et al., 2006).
6. Puyo, Pastaza — eastern-slope race records (Brown, 1981; Jiggins et al., 2008).

## 4. Range and habitat

- ***H. erato etylus*** — Eastern Andean slope and valleys of Ecuador (Napo, Pastaza, Zamora Chinchipe), ~1,000–1,400 m humid foothill/montane forest; putative hybrid race *notabilis* × *lativitta* (Jiggins et al., 2008); type locality Guadalquiza (Salvin, 1871).
- ***H. erato demophoon*** — Both coasts of Costa Rica (type locality Nicaragua) through Panama; open and disturbed habitats, forest edges, sea level–~1,500 m; postman mimicry ring with *H. melpomene rosina* in Panama (Hill et al., 2013; Kronforst & Gilbert, 2008; Mackay-Smith et al., 2026).
- ***H. himera*** — Endemic to seasonally dry forest and thorn scrub of southern Ecuador (Loja, El Oro: Vilcabamba, Catamayo, Río Puyango valleys) and northern Peru (Chinchipe/Marañón valleys), ~1,000–1,600 m; *Acacia macracantha* thorn scrub and xerophytic vegetation, strong diurnal temperature swings; larvae on *Passiflora* subg. *Decaloba*, adults on *Lantana*; parapatric contact zones with *H. erato* races track the wet↔dry forest ecotone (Borrero et al., 2026; Jiggins et al., 1996; Van Belleghem et al., 2021).
- *himera* × *etylus* proximity: humid *etylus* valleys (Nangaritza, Zamora) lie immediately adjacent to *himera* dry valleys (Catamayo/Zamora); nearest localities ≈ 42 km apart.

## 5. Figure design

- Dots only (no range polygons). Base: SRTM hillshade (2–98% quantile stretch) + 500 m contours, `elevatr`/`terra`; bathymetry masked, ocean white; no political borders.
- Single continuous panel, extent −90.7…−74.7° lon, −7.7…10.8° lat; the blank white Pacific carries the world locator inset (left-centre, `rnaturalearth` ne_110, grey rectangle = map extent), 300 km scale bar (bottom-left) and north arrow (top-left).
- Plotted symbols: filled circle = GBIF occurrence record; open square = iNaturalist research-grade observation. Literature localities (§3) are not plotted. Fill colour = taxon.
- Colours: demophoon `#994F00`, etylus `#006CD1` (himera-pair palette of `ani-analysis/plot-butterfly-ani.R`); himera `#D0D55C` (yellow, `mc` palette of the same script).

## 6. References (APA 7th)

Borrero, J., Rivas-Sánchez, D. F., Wright, D. S., Bacquet, C. N., Montgomery, S. H., Keller, A., & Merrill, R. M. (2026). Sexual dimorphism in pollen foraging and sensory traits in *Heliconius* butterflies. *Behavioral Ecology, 37*(3), arag014. https://doi.org/10.1093/beheco/arag014

Brown, K. S. (1981). The biology of *Heliconius* and related genera. *Annual Review of Entomology, 26*, 427–457. https://doi.org/10.1146/annurev.en.26.010181.002235

Global Biodiversity Information Facility. (2026). *Occurrence records of Heliconius himera (taxonKey 1900219) and Heliconius erato demophoon (taxonKey 11869372)* [Data set]. GBIF. Retrieved October 21, 2026, from https://www.gbif.org/occurrence/search?taxon_key=1900219

Hill, R., Gilbert, L., & Kronforst, M. (2013). Cryptic genetic and wing pattern diversity in a mimetic *Heliconius* butterfly. *Molecular Ecology, 22*(10), 2760–2770. https://doi.org/10.1111/mec.12290

iNaturalist. (2026). *Research-grade observations of Heliconius erato etylus (Taxon 1136782), H. erato demophoon (Taxon 1552759) and H. himera (Taxon 790498)* [Occurrence data]. iNaturalist.org. Retrieved October 21, 2026, from https://www.inaturalist.org/taxa/1136782-Heliconius-erato-etylus

Jiggins, C. D., McMillan, W. O., Neukirchen, W., & Mallet, J. (1996). What can hybrid zones tell us about speciation? The case of *Heliconius erato* and *H. himera* (Lepidoptera: Nymphalidae). *Biological Journal of the Linnean Society, 59*(3), 221–242. https://doi.org/10.1111/j.1095-8312.1996.tb01464.x

Jiggins, C. D., Salazar, C., Linares, M., & Mavarez, J. (2008). Hybrid trait speciation and *Heliconius* butterflies. *Philosophical Transactions of the Royal Society B, 363*(1508), 3047–3054. https://doi.org/10.1098/rstb.2008.0065

Kapan, D. D., Flanagan, N. S., Tobler, A., Papa, R., Reed, R. D., Acevedo Gonzalez, J., Ramirez Restrepo, M., Martinez, L., Maldonado, K., Ritschoff, C., Heckel, D. G., & McMillan, W. O. (2006). Localization of Müllerian mimicry genes on a dense linkage map of *Heliconius erato*. *Genetics, 173*(2), 735–757. https://doi.org/10.1534/genetics.106.057166

Kronforst, M. R., & Gilbert, L. E. (2008). The population genetics of mimetic diversity in *Heliconius* butterflies. *Proceedings of the Royal Society B, 275*(1634), 493–500. https://doi.org/10.1098/rspb.2007.1378

Mackay-Smith, A., Alexander, G. C., McMillan, W. O., Martin, A., & Wray, G. A. (2026). De novo assembly of five highly contiguous *Heliconius* butterfly genomes with long-read sequencing alone. *Genome Biology and Evolution, 18*(7), evag171. https://doi.org/10.1093/gbe/evag171

Salvin, O. (1871). Descriptions of new species of diurnal Lepidoptera found in Central America. *Annals and Magazine of Natural History, 4*(21), 414.

Van Belleghem, S. M., Cole, J. M., Montejo-Kovacevich, G., Bacquet, C. N., McMillan, W. O., Papa, R., & Counterman, B. A. (2021). Selection and isolation define a heterogeneous divergence landscape between hybridizing *Heliconius* butterflies. *Evolution, 75*(9), 2251–2268. https://doi.org/10.1111/evo.14272

### Web resources (no DOI)

ETH Zürich Entomological Collection. (2023). *Heliconiini collection records AZ23-0322–AZ23-0324* [Specimen labels]. https://lepidoptera.phys.ethz.ch/Collection-Heliconiini.html

Funet. (2024). *Lepidoptera and some other life forms: Heliconius* [Taxonomic index, type localities]. https://ftp.funet.fi/index/Tree_of_life/insecta/lepidoptera/ditrysia/papilionoidea/nymphalidae/heliconiinae/heliconius/

Lepiforum e.V. (2025). *Heliconius erato* [Subspecies records and photographs, LepiWiki]. https://lepiforum.org/wiki/page/Heliconius_erato

Tree of Life Web Project. (2008). *Heliconius erato* [Race distribution map]. https://tolweb.org/Heliconius_erato/72242
