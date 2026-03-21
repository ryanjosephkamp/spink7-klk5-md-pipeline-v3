# SPINK7-KLK5 MD Pipeline — Figure Generation Makefile
PYTHON := python
FIGURES_DIR := figures
ANALYSIS_DIR := data/analysis

.PHONY: figures clean-figures help

help:
	@echo "Available targets:"
	@echo "  make figures        — Run analysis and generate all publication figures"
	@echo "  make clean-figures  — Remove all generated figures"
	@echo "  make help           — Show this help message"

figures: $(FIGURES_DIR)/fig3_protein_complex.png \
         $(FIGURES_DIR)/fig4_pmf_profile.png \
         $(FIGURES_DIR)/fig5_simulation_timeseries.png
	@echo "All publication figures generated in $(FIGURES_DIR)/"

$(FIGURES_DIR)/fig3_protein_complex.png $(FIGURES_DIR)/fig4_pmf_profile.png $(FIGURES_DIR)/fig5_simulation_timeseries.png: scripts/generate_figures.py
	$(PYTHON) scripts/generate_figures.py

clean-figures:
	rm -f $(FIGURES_DIR)/fig3_protein_complex.{png,svg,pdf}
	rm -f $(FIGURES_DIR)/fig4_pmf_profile.{png,svg,pdf}
	rm -f $(FIGURES_DIR)/fig5_simulation_timeseries.{png,svg,pdf}
