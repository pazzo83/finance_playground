using Base.Dates
using Requests
using DataFrames
using TimeSeries

abstract DataSource
immutable YahooDataSource <: DataSource end

# type YahooFinanceQuery{S <: AbstractString}
#   query::S
#   responseJSON::S
# end
#
# type YahooDataSet{S <: AbstractString}
#   rawJson::S
#   dates::Vector{Date}
#   lows::Vector{Float64}
#   highs::Vector{Float64}
#   opens::Vector{Float64}
#   closes::Vector{Float64}
#   volumes::Vector{Float64}
#   adjCloses::Vector{Float64}
#   symbs::Vector{S}
# end

# function get_historical_data(::YahooDataSource, startDate::Date, endDate::Date, symbol::AbstractString)
#   # static dates right now
#   # date1 = Date(2009, 9, 11)
#   # date2 = Date(2010, 3, 10)
#   # symb = "YHOO"
#
#   # construct query
#   url = "https://query.yahooapis.com/v1/public/yql?q="
#
#   query = "select%20*%20from%20yahoo.finance.historicaldata%20where%20symbol%20%3D%20%22$symbol%22%20and%20startDate%20%3D%20%22$(string(startDate))%22%20and%20endDate%20%3D%20%22$(string(endDate))%22"
#
#   resp = get(string(url, query, "&format=json&diagnostics=true&env=store%3A%2F%2Fdatatables.org%2Falltableswithkeys&callback="))
#
#   respJson = Requests.json(resp)["query"]["results"]["quote"]
#
#   n = length(respJson)
#   dates = Vector{Date}(n)
#   lows = Vector{Float64}(n)
#   highs = Vector{Float64}(n)
#   opens = Vector{Float64}(n)
#   closes = Vector{Float64}(n)
#   volumes = Vector{Float64}(n)
#   adjCloses = Vector{Float64}(n)
#   symbs = Vector{ASCIIString}(n)
#
#   for i in eachindex(respJson)
#     dates[i] = Date(respJson[i]["Date"])
#     lows[i] = float(respJson[i]["Low"])
#     highs[i] = float(respJson[i]["High"])
#     opens[i] = float(respJson[i]["Open"])
#     closes[i] = float(respJson[i]["Close"])
#     volumes[i] = float(respJson[i]["Volume"])
#     adjCloses[i] = float(respJson[i]["Adj_Close"])
#     symbs[i] = respJson[i]["Symbol"]
#   end
#
#   df = DataFrame(Date = dates, Low = lows, High = highs, Close = closes, Volume = volumes, AdjClose = adjCloses, Symbol = symbs)
#
#   return df
# end

function get_historical_data(::YahooDataSource, startDate::Date, endDate::Date, symbol::AbstractString)
  startDate < endDate || error("dates are wrong")

  url = "http://ichart.finance.yahoo.com/table.csv?s="

  query = "$symbol&d=$(Int(Month(endDate)) - 1)&e=$(Int(Day(endDate)))&f=$(Int(Year(endDate)))&g=d&a=$(Int(Month(startDate)) - 1)&b=$(Int(Day(startDate)))&c=$(Int(Year(startDate)))"

  # resp = get(url, query, "&ignore.csv")
  #
  # respCSV = readall(resp)

  urlFull = string(url, query, "&ignore.csv")

  # println(urlFull)

  r = get(urlFull)

  # println(r.response)

  # df = readtable(get_streaming(string(url, query, "&ignore.csv")))

  df = readtable(IOBuffer(readall(r)))
  df[:Date] = Date(df[:Date])

  return df
end

function read_dax_data()
  DAX = get_historical_data(YahooDataSource(), Date(2004, 9, 30), Date(2014, 9, 30), "^GDAXI")

  # convert to time series for calc
  ts = TimeArray(DAX[:Date].data, DAX[:Adj_Close], ["returns"])
  tsCalc = log(ts["returns"] ./ lag(ts["returns"], padding=true))
  DAX[:returns] = tsCalc.values
  DAX[isnan(DAX[:returns]), :returns] = 0.0

  # Realized volatility (e.g. as defined for variance swaps)
  DAX[:rea_var] = 252 * cumsum(DAX[:returns] .^ 2) ./ range(1, size(DAX)[1])
  DAX[:rea_vol] = sqrt(DAX[:rea_var])

  return DAX
end
