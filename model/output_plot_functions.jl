module output_plot_functions

using Plots
using Statistics


function plot_dispatch(MODEL_INPUT, RESULTS_ALL)

    
end


function plot_lmp_aggregated(MODEL_INPUT, RESULTS_ALL)

    cols_prices = ["Max. LMP", "95% Percentile", "75%-Quantile", "Median LMP", "25%-Quantile", "5% Percentile", "Min. LMP"]
    
    size(RESULTS_ALL.PRICES)
    data_sorted = sort(RESULTS_ALL.PRICES, dims=1)
    relevant_rows_zero = round.(Int, size(RESULTS_ALL.PRICES, 1)*[1, 0.95, 0.75, 0.5, 0.25, 0.05, 0.0])
    relevant_rows = relevant_rows_zero
    relevant_rows[end] = 1

    data_vals = data_sorted[relevant_rows,:]'

    bla = plot(data_vals, linetype=:steppre, labels=permutedims(cols_prices),legend=:best, title="LMP", ylabel="€/MWh", xlabel="hour")
    display(bla)

end

end

