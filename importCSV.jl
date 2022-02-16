using CSV
using DataFrames
using Statistics



function importCSV()
    df = CSV.read("C:\\Users\\Jake\\Downloads\\Discovering-Multipoles-master\\Discovering-Multipoles-master\\New\\multipolesCSV.csv", DataFrame,header = false)
    M = Matrix(df)
    println(size(M))
    C = cov(M)
    return C
end
