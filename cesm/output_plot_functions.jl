module output_plot_functions

using Plots
using Statistics
using StatsPlots


function plot_generation_stack(MODEL_INPUT, RESULTS_ALL)

    
    data_matrix = vcat(
        sum(RESULTS_ALL.RENEWABLE_FEED_IN,dims=1),    
        sum(RESULTS_ALL.GEN_EL_ONLY,dims=1),
        sum(RESULTS_ALL.GEN_EL_CHP,dims=1),
        sum(RESULTS_ALL.GEN_HY_ROR,dims=1),
        sum(RESULTS_ALL.GEN_EL_BIO,dims=1),
    )


    figure1 = groupedbar(data_matrix'./1e3,
        bar_position = :stack,
        bar_width=1,
        color=permutedims(["green", "orange", "purple", "blue", "brown"]),
        label=permutedims(["Wind & PV", "Gas", "CHP", "HY_ROR", "Biomass"]),
        xlabel="Timesteps",
        ylabel="Power in GW",
        legend=:best
        )
    display(figure1)



end


function plot_lmp_aggregated(MODEL_INPUT, RESULTS_ALL)

    cols_prices = ["Max. LMP", "95% Percentile", "75%-Quantile", "Median LMP", "25%-Quantile", "5% Percentile", "Min. LMP"]
    
    size(RESULTS_ALL.PRICES)
    data_sorted = sort(RESULTS_ALL.PRICES, dims=1)
    relevant_rows_zero = round.(Int, size(RESULTS_ALL.PRICES, 1)*[1, 0.95, 0.75, 0.5, 0.25, 0.05, 0.0])
    relevant_rows = relevant_rows_zero
    relevant_rows[end] = 1

    data_vals = data_sorted[relevant_rows,:]'

    figure2 = plot(data_vals, linetype=:steppre, legend=:best, labels=permutedims(cols_prices), 
        title="Locational Marginal Prices", ylabel="Price in €/MWh", xlabel="Timesteps")
    display(figure2)

end


end

