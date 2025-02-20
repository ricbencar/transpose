// --------------------------------------------------------------------
// OFFSHORE-TO-NEARSHORE WAVE TRANSFORMATION
// (Simplified Approach Using Linear Wave Theory and Hybrid Circular Statistics)
//
// Overview:
//   This program processes wave data from an input CSV file and computes
//   nearshore wave parameters at a specified depth. It generates two output files:
//     - "output.csv" – Contains the computed nearshore wave parameters.
//     - "report.txt" – Provides descriptive statistics for both the input and
//                      computed variables.
//   The report includes:
//     * The command line used to invoke the program.
//     * Descriptive statistics for each variable (count, mean, standard deviation,
//       minimum, maximum, median, and percentiles at 1%, 10%, 25%, 50%, 75%, 90%, 
//       and 99%).
//     * A table of annual maxima for swh_offshore and swh_local, with a final row
//       indicating the overall maximum values.
//
//   For directional wave data (mwd_offshore and mwd_local), a hybrid approach is used:
//     - The circular mean and circular standard deviation are computed using the
//       unit-vector method.
//     - The minimum, maximum, median, and quantiles are calculated using ordinary
//       linear statistics on the wrapped angles (in [0,360)).
//
// USAGE:
//   ./transpose input_csv coast_dir depth_d
//
//   Where:
//     input_csv : CSV input file containing at least the following columns:
//                 datetime, swh, mwd, pp1d (additional columns are ignored)
//     coast_dir : Coastline orientation in degrees (clockwise from North)
//     depth_d   : Local depth in meters
//
// EXPECTED CSV INPUT FORMAT (comma-separated):
//     datetime, swh, mwd, pp1d, [additional columns ignored]
//
// OUTPUT CSV FORMAT (comma-separated):
//     datetime,swh_offshore,mwd_offshore,pp1d,L0,L,kh,alpha_offshore,
//     alpha_local,swh_local,mwd_local,Ks,Kr,Hb
//
// Explanation of computed parameters:
//     L0             : Deep-water wavelength, calculated as (g * T²) / (2π)
//     L              : Local wavelength, solved via Newton-Raphson from
//                      L = L0 * tanh((2π * depth_d) / L)
//     kh             : Product of the wave number (k = 2π / L) and local depth (h)
//     alpha_offshore : Offshore wave approach angle relative to the coastline
//     alpha_local    : Local wave angle after refraction
//     mwd_local      : Local mean wave direction, adjusted from the offshore mwd
//     Ks             : Shoaling coefficient
//     Kr             : Refraction coefficient
//     Hb             : Breaking wave height (per Miche, 1944), computed as
//                      Hb = 0.142 * L * tanh((2π * depth_d) / L)
//     swh_local      : Local significant wave height, computed as the minimum of
//                      (swh * Ks * Kr) and Hb
//
//   Note: Waves arriving from directions between coast_dir and coast_dir+180°
//         (i.e., from the land side) are set to zero.
//
// Report File Details:
//   The report.txt file includes:
//     - The exact command line used to run the program.
//     - Detailed descriptive statistics for each variable, including count,
//       mean, standard deviation, minimum, maximum, median, and percentiles at
//       1%, 10%, 25%, 50%, 75%, 90%, and 99%.
//     - A table displaying the annual maxima for swh_offshore and swh_local,
//       with the final row indicating the overall maximum values.
//
// Compilation Details:
//   To compile the program, use the following command:
//
//     g++ -O3 -fopenmp -march=native -std=c++17 -Wall -Wextra -pedantic
//         -Wconversion -Wsign-conversion -static -static-libgcc -static-libstdc++
//         -o transpose transpose.cpp
//
//   Explanation of compile options:
//     - -O3                   : Enables high-level optimizations for maximum performance.
//     - -fopenmp              : Enables OpenMP support for multi-threading.
//     - -march=native         : Optimizes the code for the architecture of the compiling machine.
//     - -std=c++17            : Uses the C++17 standard.
//     - -Wall -Wextra -pedantic: Activates a broad set of compiler warnings to ensure code quality.
//     - -Wconversion          : Warns about implicit type conversions.
//     - -Wsign-conversion     : Warns about implicit sign conversions.
//     - -static, -static-libgcc, -static-libstdc++: Links libraries statically, enhancing portability.
// --------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>   // For std::setprecision, std::setw, std::left
#include <limits>
#include <cstdlib>
#include <map>

using namespace std;

// --------------------------------------------------------------------
// Physical constants and iteration parameters
// --------------------------------------------------------------------
static const long double G = 9.81L;                  // Acceleration due to gravity (m/s²)
static const long double PI = 3.141592653589793238L; // Mathematical constant π
static const int MAX_ITER = 20;                      // Maximum iterations for Newton-Raphson solver
static const long double TOLERANCE = 1e-12L;         // Convergence tolerance

// --------------------------------------------------------------------
// Utility functions for angle conversion
// --------------------------------------------------------------------
inline long double deg2rad(long double deg)
{
    return deg * (PI / 180.0L);
}

inline long double rad2deg(long double rad)
{
    return rad * (180.0L / PI);
}

// --------------------------------------------------------------------
// Calculate offshore obliquity (alpha_offshore in degrees).
// --------------------------------------------------------------------
long double calcalpha_offshore_deg(long double mwd_deg, long double coast_deg)
{
    long double crest = mwd_deg - 90.0L;
    // Normalize crest to [0, 360)
    crest = fmodl(crest + 360.0L, 360.0L);
    long double diff = fabsl(crest - coast_deg);
    diff = fmodl(diff, 360.0L);
    if (diff > 180.0L)
        diff = 360.0L - diff;
    if (diff > 90.0L)
        diff = 180.0L - diff;
    return diff;
}

// --------------------------------------------------------------------
// Deep-water wavelength calculation.
// --------------------------------------------------------------------
long double deepWaterLength(long double T)
{
    return (G * T * T) / (2.0L * PI);
}

// --------------------------------------------------------------------
// Local wavelength computation using Newton-Raphson.
// --------------------------------------------------------------------
long double localWavelength(long double T, long double depth)
{
    if (T <= 0.0L || depth <= 0.0L)
        return 0.0L;

    const long double L0 = deepWaterLength(T);
    long double k0 = (2.0L * PI * depth) / L0;
    // Use a simple initial guess based on deep-water wavelength
    long double L = L0 * tanhl(sinhl(sqrt(k0))) / tanhl(T);

    for (int i = 0; i < MAX_ITER; ++i)
    {
        long double x = (2.0L * PI * depth) / L;
        long double th = tanhl(x);
        long double F = L - L0 * th;
        long double sech2 = 1.0L - th * th;
        long double dF = 1.0L + L0 * (2.0L * PI * depth) / (L * L) * sech2;
        if (fabsl(dF) < TOLERANCE)
            break;
        long double L_new = L - F / dF;
        if (fabsl(L_new - L) < TOLERANCE)
        {
            L = L_new;
            break;
        }
        L = L_new;
    }
    return (L < 0.0L ? 0.0L : L);
}

// --------------------------------------------------------------------
// Shoaling coefficient computation.
// --------------------------------------------------------------------
long double shoalingCoefficient(long double k, long double depth)
{
    if (k <= 0.0L || depth <= 0.0L)
        return 1.0L;
    long double kh = k * depth;
    long double sinh_kh = sinhl(kh);
    long double cosh_kh = coshl(kh);
    long double denom = kh + sinh_kh * cosh_kh;
    return (denom < TOLERANCE ? 1.0L : cosh_kh / sqrtl(denom));
}

// --------------------------------------------------------------------
// Structure for descriptive statistics
// --------------------------------------------------------------------
struct DescriptiveStats
{
    size_t count;
    long double mean;    // For directional data, circular mean will be used
    long double stddev;  // For directional data, circular std. dev. will be used
    long double min;
    long double p1;     // 1st percentile
    long double p10;    // 10th percentile
    long double p25;    // 25th percentile
    long double median; // 50th percentile
    long double p75;    // 75th percentile
    long double p90;    // 90th percentile
    long double p99;    // 99th percentile
    long double max;
};

// --------------------------------------------------------------------
// Compute descriptive statistics (for linear data).
// --------------------------------------------------------------------
DescriptiveStats computeStats(const vector<long double> &data)
{
    DescriptiveStats stats = {0, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
    if (data.empty())
        return stats;
    stats.count = data.size();
    long double sum = 0.0L;
    stats.min = data[0];
    stats.max = data[0];
    for (long double x : data)
    {
        sum += x;
        if (x < stats.min)
            stats.min = x;
        if (x > stats.max)
            stats.max = x;
    }
    stats.mean = sum / stats.count;
    long double variance = 0.0L;
    for (long double x : data)
    {
        long double diff = x - stats.mean;
        variance += diff * diff;
    }
    stats.stddev = (stats.count > 1) ? sqrtl(variance / (stats.count - 1)) : 0.0L;

    vector<long double> sortedData = data;
    sort(sortedData.begin(), sortedData.end());
    auto getPercentile = [&](long double p) -> long double {
        long double pos = p * (sortedData.size() - 1);
        size_t idx = (size_t)floorl(pos);
        long double frac = pos - idx;
        return (idx + 1 < sortedData.size())
            ? sortedData[idx] * (1.0L - frac) + sortedData[idx + 1] * frac
            : sortedData[idx];
    };
    stats.p1     = getPercentile(0.01L);
    stats.p10    = getPercentile(0.10L);
    stats.p25    = getPercentile(0.25L);
    stats.median = getPercentile(0.50L);
    stats.p75    = getPercentile(0.75L);
    stats.p90    = getPercentile(0.90L);
    stats.p99    = getPercentile(0.99L);
    return stats;
}

// --------------------------------------------------------------------
// Compute hybrid circular statistics for directional data.
// For mwd_offshore and mwd_local, the unit-vector method is used
// to compute the circular mean and standard deviation, ordinary (linear)
// statistics are used on the "wrapped" angles (in [0,360)),
// for min, max, median, and quantiles.
// --------------------------------------------------------------------
DescriptiveStats computeHybridCircularStats(const vector<long double> &data)
{
    // Wrap each angle into [0,360)
    vector<long double> wrapped;
    wrapped.reserve(data.size());
    for (long double d : data) {
        long double w = fmodl(d, 360.0L);
        if (w < 0)
            w += 360.0L;
        wrapped.push_back(w);
    }
    
    // Compute linear statistics on the wrapped values
    DescriptiveStats linearStats = computeStats(wrapped);
    
    // Compute circular mean and circular standard deviation via unit-vector method.
    long double sumSin = 0.0L, sumCos = 0.0L;
    for (long double d : wrapped) {
        long double rad = deg2rad(d);
        sumSin += sinl(rad);
        sumCos += cosl(rad);
    }
    long double meanRad = atan2l(sumSin, sumCos);
    if (meanRad < 0)
        meanRad += 2.0L * PI;
    long double circMean = rad2deg(meanRad);
    long double R = sqrtl((sumCos / wrapped.size()) * (sumCos / wrapped.size()) +
                          (sumSin / wrapped.size()) * (sumSin / wrapped.size()));
    long double circStd = (R < TOLERANCE) ? 180.0L : rad2deg(sqrtl(-2.0L * logl(R)));
    
    // Build the hybrid statistics: use circular mean and std. dev. and linear quantiles.
    DescriptiveStats hybrid = linearStats;
    hybrid.mean = circMean;
    hybrid.stddev = circStd;
    return hybrid;
}

// --------------------------------------------------------------------
// Main function
// --------------------------------------------------------------------
int main(int argc, char *argv[])
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc != 4)
    {
        cerr << "Usage: " << argv[0] << " input_csv coast_dir depth_d\n"
             << "Example: " << argv[0] << " input.csv 10 5\n";
        return 1;
    }
    string input_csv = argv[1];
    long double coast_dir, depth_d;
    try
    {
        coast_dir = stold(argv[2]);
        depth_d = stold(argv[3]);
    }
    catch (...)
    {
        cerr << "Error: cannot parse coast_dir or depth_d.\n";
        return 1;
    }
    if (depth_d <= 0.0L)
    {
        cerr << "Invalid depth value.\n";
        return 1;
    }

    ifstream inFile(input_csv);
    if (!inFile.is_open())
    {
        cerr << "ERROR: unable to open input file " << input_csv << "\n";
        return 1;
    }
    ofstream outFile("output.csv");
    if (!outFile.is_open())
    {
        cerr << "ERROR: unable to create output.csv\n";
        return 1;
    }
    // Write output CSV header
    outFile << "datetime,"
            << "swh_offshore,mwd_offshore,pp1d,"
            << "L0,L,kh,alpha_offshore,alpha_local,"
            << "swh_local,mwd_local,Ks,Kr,Hb\n";

    // Prepare arrays for statistical data (13 columns)
    const size_t NUM_COLS = 13;
    vector<vector<long double>> statsData(NUM_COLS);
    vector<string> varNames = {
        "swh_offshore",   // 0
        "mwd_offshore",   // 1
        "pp1d",           // 2
        "L0",             // 3
        "L",              // 4
        "kh",             // 5
        "alpha_offshore", // 6
        "alpha_local",    // 7
        "swh_local",      // 8
        "mwd_local",      // 9
        "Ks",             // 10
        "Kr",             // 11
        "Hb"              // 12
    };
    for (auto &col : statsData)
        col.reserve(1000000);

    // Maps to record annual maximum values for swh_offshore and swh_local.
    map<string, long double> annualMaxOffshore;
    map<string, long double> annualMaxLocal;

    // Read and discard CSV header
    string header;
    getline(inFile, header);

    // Read remaining lines into a vector
    vector<string> lines;
    string line;
    while (getline(inFile, line))
    {
        if (!line.empty())
            lines.push_back(line);
    }
    inFile.close();

    // Sort lines by datetime (first token) and remove duplicate datetime entries.
    sort(lines.begin(), lines.end(), [](const string &a, const string &b) {
        size_t posA = a.find(',');
        size_t posB = b.find(',');
        string dtA = (posA == string::npos) ? a : a.substr(0, posA);
        string dtB = (posB == string::npos) ? b : b.substr(0, posB);
        return dtA < dtB;
    });
    vector<string> uniqueLines;
    for (const auto &l : lines)
    {
        size_t pos = l.find(',');
        string dt = (pos == string::npos) ? l : l.substr(0, pos);
        if (!uniqueLines.empty())
        {
            size_t posLast = uniqueLines.back().find(',');
            string lastDt = (posLast == string::npos) ? uniqueLines.back() : uniqueLines.back().substr(0, posLast);
            if (lastDt == dt)
                continue; // duplicate found
        }
        uniqueLines.push_back(l);
    }

    // Process each unique line
    for (const auto &l : uniqueLines)
    {
        size_t pos = 0, next;
        // Parse datetime
        next = l.find(',', pos);
        if (next == string::npos)
            continue;
        string datetime = l.substr(pos, next - pos);
        pos = next + 1;
        // Parse swh (this is swh_offshore)
        next = l.find(',', pos);
        if (next == string::npos)
            continue;
        string sSwh = l.substr(pos, next - pos);
        pos = next + 1;
        // Parse mwd
        next = l.find(',', pos);
        if (next == string::npos)
            continue;
        string sMwd = l.substr(pos, next - pos);
        pos = next + 1;
        // Parse pp1d
        next = l.find(',', pos);
        string sPp1d = (next == string::npos) ? l.substr(pos) : l.substr(pos, next - pos);

        long double swh, mwd, pp1d;
        try
        {
            swh = stold(sSwh);
            mwd = stold(sMwd);
            pp1d = stold(sPp1d);
        }
        catch (...)
        {
            continue; // Skip invalid numeric lines
        }
        
        // We update annual maximum only if datetime is at least 4 characters.
        bool updateMaps = (datetime.size() >= 4);

        // Branch 1: Invalid numeric values
        if (swh <= 0.0L || pp1d <= 0.0L)
        {
            outFile << datetime << "," << swh << "," << mwd << "," << pp1d;
            for (int i = 0; i < 9; i++)
                outFile << ",NaN";
            outFile << "\n";
            vector<long double> record = {swh, mwd, pp1d, 0.0L, 0.0L, 0.0L, 0.0L,
                                           0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
            for (size_t i = 0; i < NUM_COLS; i++)
                statsData[i].push_back(record[i]);
            continue;
        }

        // Compute relative direction
        long double relativeDir = fmodl((mwd - coast_dir) + 360.0L, 360.0L);

        // Branch 2: Waves arriving from the land side (relativeDir < 180)
        if (relativeDir < 180.0L)
        {
            // swh_offshore equals swh; swh_local is set to 0.0 in this branch.
            long double swh_local = 0.0L;
            outFile << datetime << "," << swh << "," << mwd << "," << pp1d;
            for (int i = 0; i < 9; i++)
                outFile << ",0.0";
            outFile << "\n";
            vector<long double> record = {swh, mwd, pp1d, 0.0L, 0.0L, 0.0L, 0.0L,
                                           0.0L, swh_local, 0.0L, 0.0L, 0.0L, 0.0L};
            for (size_t i = 0; i < NUM_COLS; i++)
                statsData[i].push_back(record[i]);
            if (updateMaps)
            {
                string year = datetime.substr(0, 4);
                // Update offshore maximum
                if (annualMaxOffshore.find(year) == annualMaxOffshore.end())
                    annualMaxOffshore[year] = swh;
                else
                    annualMaxOffshore[year] = max(annualMaxOffshore[year], swh);
                // Update local maximum (0.0 in this branch)
                if (annualMaxLocal.find(year) == annualMaxLocal.end())
                    annualMaxLocal[year] = swh_local;
                else
                    annualMaxLocal[year] = max(annualMaxLocal[year], swh_local);
            }
            continue;
        }

        // Branch 3: Valid offshore wave case
        long double L0 = deepWaterLength(pp1d);
        long double L = localWavelength(pp1d, depth_d);
        long double kLocal = (L > 1e-12L) ? (2.0L * PI / L) : 0.0L;
        long double kh = kLocal * depth_d;
        long double alpha_offshore_deg = calcalpha_offshore_deg(mwd, coast_dir);
        long double alpha_offshore_rad = deg2rad(alpha_offshore_deg);
        long double sinAlpha = sinl(alpha_offshore_rad) * tanhl(kh);
        if (fabsl(sinAlpha) > 1.0L)
            sinAlpha = (sinAlpha > 0.0L ? 1.0L : -1.0L);
        long double alpha_local_rad = asinl(sinAlpha);
        long double alpha_local_deg = rad2deg(alpha_local_rad);
        long double mwd_local = fmodl(mwd - (alpha_offshore_deg - alpha_local_deg) + 360.0L, 360.0L);
        long double Ks = shoalingCoefficient(kLocal, depth_d);
        long double cosalpha_offshore = cosl(alpha_offshore_rad);
        long double cosAlphaLocal = cosl(alpha_local_rad);
        long double Kr = (fabsl(cosAlphaLocal) > 1e-12L) ? sqrtl(cosalpha_offshore / cosAlphaLocal) : 1.0L;
        long double Hb = (L > 1e-12L) ? 0.142L * L * tanhl((2.0L * PI * depth_d) / L) : 0.0L;
        long double swh_local = swh * Ks * Kr;
        if (Hb > 0.0L && swh_local > Hb)
            swh_local = Hb;
        outFile << datetime << "," << swh << "," << mwd << "," << pp1d << ","
                << L0 << "," << L << "," << kh << ","
                << alpha_offshore_deg << "," << alpha_local_deg << ","
                << swh_local << "," << mwd_local << ","
                << Ks << "," << Kr << "," << Hb << "\n";
        vector<long double> record = {swh, mwd, pp1d, L0, L, kh,
                                      alpha_offshore_deg, alpha_local_deg, swh_local,
                                      mwd_local, Ks, Kr, Hb};
        for (size_t i = 0; i < NUM_COLS; i++)
            statsData[i].push_back(record[i]);
        if (updateMaps)
        {
            string year = datetime.substr(0, 4);
            // Update offshore maximum (swh_offshore equals swh)
            if (annualMaxOffshore.find(year) == annualMaxOffshore.end())
                annualMaxOffshore[year] = swh;
            else
                annualMaxOffshore[year] = max(annualMaxOffshore[year], swh);
            // Update local maximum
            if (annualMaxLocal.find(year) == annualMaxLocal.end())
                annualMaxLocal[year] = swh_local;
            else
                annualMaxLocal[year] = max(annualMaxLocal[year], swh_local);
        }
    }
    outFile.close();

    // ----------------------------------------------------------------
    // Create report.txt with descriptive statistics and the annual
    // maxima table (for swh_offshore and swh_local) with overall maxima.
    // ----------------------------------------------------------------
    const int LINE_WIDTH = 100;                    // Total characters per line
    const int VARIABLES_PER_ROW = 3;               // Display 3 variables side by side
    int colWidth = LINE_WIDTH / VARIABLES_PER_ROW; // Width per variable

    ofstream reportFile("report.txt");
    if (!reportFile.is_open())
    {
        cerr << "ERROR: unable to create report.txt\n";
        return 1;
    }
    // Write the command line used at the start of the report file
    reportFile << "Command line: " << argv[0] << " " << argv[1] << " " << argv[2] << " " << argv[3] << "\n\n";
    reportFile << "Descriptive Statistics Report\n";
    reportFile << string(LINE_WIDTH, '=') << "\n\n";

    // Define the list of statistic labels for each variable.
    vector<string> statLabels = {"Count", "Mean", "Std. Dev.", "Min", "P1", "P10", "P25", "Median", "P75", "P90", "P99", "Max"};
    for (size_t i = 0; i < NUM_COLS; i += VARIABLES_PER_ROW)
    {
        // Print variable names for this group
        for (size_t j = 0; j < VARIABLES_PER_ROW && (i + j) < NUM_COLS; ++j)
            reportFile << left << setw(colWidth) << varNames[i + j];
        reportFile << "\n";
        vector<DescriptiveStats> groupStats;
        for (size_t j = 0; j < VARIABLES_PER_ROW && (i + j) < NUM_COLS; ++j) {
            // For directional variables, use the hybrid method.
            if (varNames[i + j] == "mwd_offshore" || varNames[i + j] == "mwd_local")
                groupStats.push_back(computeHybridCircularStats(statsData[i + j]));
            else
                groupStats.push_back(computeStats(statsData[i + j]));
        }
        // For each statistic, print the value for each variable
        for (const auto &label : statLabels)
        {
            for (size_t j = 0; j < groupStats.size(); ++j)
            {
                ostringstream oss;
                if (label == "Count")
                    oss << label << ": " << groupStats[j].count;
                else if (label == "Mean")
                    oss << label << ": " << fixed << setprecision(10) << groupStats[j].mean;
                else if (label == "Std. Dev.")
                    oss << label << ": " << fixed << setprecision(10) << groupStats[j].stddev;
                else if (label == "Min")
                    oss << label << ": " << fixed << setprecision(10) << groupStats[j].min;
                else if (label == "P1")
                    oss << label << ": " << fixed << setprecision(10) << groupStats[j].p1;
                else if (label == "P10")
                    oss << label << ": " << fixed << setprecision(10) << groupStats[j].p10;
                else if (label == "P25")
                    oss << label << ": " << fixed << setprecision(10) << groupStats[j].p25;
                else if (label == "Median")
                    oss << label << ": " << fixed << setprecision(10) << groupStats[j].median;
                else if (label == "P75")
                    oss << label << ": " << fixed << setprecision(10) << groupStats[j].p75;
                else if (label == "P90")
                    oss << label << ": " << fixed << setprecision(10) << groupStats[j].p90;
                else if (label == "P99")
                    oss << label << ": " << fixed << setprecision(10) << groupStats[j].p99;
                else if (label == "Max")
                    oss << label << ": " << fixed << setprecision(10) << groupStats[j].max;
                reportFile << left << setw(colWidth) << oss.str();
            }
            reportFile << "\n";
        }
        reportFile << string(LINE_WIDTH, '-') << "\n";
    }
    
    // Append the Annual Maxima table for swh_offshore and swh_local.
    reportFile << "\nAnnual Maxima of swh_offshore and swh_local:\n";
    // Table header
    int tableColWidth = 30;
    reportFile << left << setw(tableColWidth) << "Year" 
               << left << setw(tableColWidth) << "swh_offshore"
               << left << setw(tableColWidth) << "swh_local" << "\n";
    long double overallOffshore = 0.0L, overallLocal = 0.0L;
    for (const auto &entry : annualMaxOffshore)
    {
        string year = entry.first;
        long double offVal = entry.second;
        long double localVal = 0.0L;
        if (annualMaxLocal.find(year) != annualMaxLocal.end())
            localVal = annualMaxLocal[year];
        reportFile << left << setw(tableColWidth) << year
                   << left << setw(tableColWidth) << offVal
                   << left << setw(tableColWidth) << localVal << "\n";
        overallOffshore = max(overallOffshore, offVal);
        overallLocal = max(overallLocal, localVal);
    }
    // Overall maxima line
    reportFile << left << setw(tableColWidth) << "Overall Max"
               << left << setw(tableColWidth) << overallOffshore
               << left << setw(tableColWidth) << overallLocal << "\n";
    
    reportFile.close();
    return 0;
}

/*
MIT License

Copyright (c) 2025 Author

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is furnished
to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/