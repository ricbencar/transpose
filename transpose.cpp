// --------------------------------------------------------------------
// OFFSHORE-TO-NEARSHORE WAVE TRANSFORMATION
// (Simplified Approach Using Linear Wave Theory)
//
// This program processes wave data from an input CSV file, computes
// nearshore wave parameters at a specified depth, and generates:
//   - "output.csv" with the computed results
//   - "report.txt" with descriptive statistics of both input and
//     computed variables. The report is formatted so that three
//     variables are shown side by side to make use of the full page width.
//
// USAGE:
//    ./transpose input_csv coast_dir depth_d
// where:
//    input_csv : CSV input file (with columns: datetime, swh, mwd, pp1d)
//    coast_dir : Coastline orientation in degrees (clockwise from North)
//    depth_d   : Local depth (meters)
//
// EXPECTED CSV INPUT FORMAT (comma-separated):
//    datetime, swh, mwd, pp1d, [additional columns ignored]
//
// OUTPUT CSV FORMAT (comma-separated):
//    datetime,swh_offshore,mwd_offshore,pp1d,L0,depth_d,L,kh,
//    alpha_offshore,alpha_local,swh_local,mwd_local,Ks,Kr,Hb
//
// Explanation of computed parameters:
//    L0             : Deep-water wavelength = g * T² / (2π)
//    depth_d        : Local depth (input parameter echoed in output)
//    L              : Local wavelength, solved from
//                     L = L0 * tanh((2π * depth_d) / L)
//    kh             : Wave number (k) times local depth (h), k = 2π / L
//    alpha_offshore : Offshore wave approach angle relative to coastline
//    alpha_local    : Local wave angle after refraction
//    mwd_local      : Local mean wave direction, adjusted from offshore mwd
//    Ks             : Shoaling coefficient
//    Kr             : Refraction coefficient
//    Hb             : Breaking wave height (Miche, 1944):
//                     Hb = 0.142 * L * tanh((2π * depth_d) / L)
//    swh_local      : Local significant wave height (min of swh * Ks * Kr and Hb)
//
// Waves arriving from directions between coast_dir and coast_dir+180°
// (i.e., from the land side) are set to zero.
//
// COMPILATION:
//    g++ -O3 -std=c++17 -Wall -Wextra -pedantic -Wconversion
//        -Wsign-conversion -o transpose transpose.cpp
// --------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip> // For std::setprecision, std::setw, std::left
#include <limits>
#include <cstdlib>

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
    long double mean;
    long double stddev;
    long double min;
    long double Q1;
    long double median;
    long double Q3;
    long double max;
};

// --------------------------------------------------------------------
// Compute descriptive statistics.
// --------------------------------------------------------------------
DescriptiveStats computeStats(const vector<long double> &data)
{
    DescriptiveStats stats = {0, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
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
    auto getPercentile = [&](long double p) -> long double
    {
        long double pos = p * (sortedData.size() - 1);
        size_t idx = (size_t)floorl(pos);
        long double frac = pos - idx;
        return (idx + 1 < sortedData.size()) ? sortedData[idx] * (1.0L - frac) + sortedData[idx + 1] * frac
                                             : sortedData[idx];
    };
    stats.median = getPercentile(0.5L);
    stats.Q1 = getPercentile(0.25L);
    stats.Q3 = getPercentile(0.75L);
    return stats;
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
    outFile << "datetime,"
            << "swh_offshore,mwd_offshore,pp1d,"
            << "L0,depth_d,L,kh,alpha_offshore,alpha_local,"
            << "swh_local,mwd_local,Ks,Kr,Hb\n";

    // Prepare arrays for statistical data
    const size_t NUM_COLS = 14;
    vector<vector<long double>> statsData(NUM_COLS);
    vector<string> varNames = {
        "swh_offshore",   // 0
        "mwd_offshore",   // 1
        "pp1d",           // 2
        "L0",             // 3
        "depth_d",        // 4
        "L",              // 5
        "kh",             // 6
        "alpha_offshore", // 7
        "alpha_local",    // 8
        "swh_local",      // 9
        "mwd_local",      // 10
        "Ks",             // 11
        "Kr",             // 12
        "Hb"              // 13
    };
    for (auto &col : statsData)
        col.reserve(1000000);

    // Discard header of input CSV
    string header;
    getline(inFile, header);

    // Use manual parsing via string::find (faster than istringstream)
    string line;
    while (getline(inFile, line))
    {
        if (line.empty())
            continue;
        size_t pos = 0, next;
        // Parse datetime
        next = line.find(',', pos);
        if (next == string::npos)
            continue;
        string datetime = line.substr(pos, next - pos);
        pos = next + 1;
        // Parse swh
        next = line.find(',', pos);
        if (next == string::npos)
            continue;
        string sSwh = line.substr(pos, next - pos);
        pos = next + 1;
        // Parse mwd
        next = line.find(',', pos);
        if (next == string::npos)
            continue;
        string sMwd = line.substr(pos, next - pos);
        pos = next + 1;
        // Parse pp1d
        next = line.find(',', pos);
        string sPp1d = (next == string::npos) ? line.substr(pos) : line.substr(pos, next - pos);

        long double swh, mwd, pp1d;
        try
        {
            swh = stold(sSwh);
            mwd = stold(sMwd);
            pp1d = stold(sPp1d);
        }
        catch (...)
        {
            continue; // Skip lines with invalid numbers
        }
        if (swh <= 0.0L || pp1d <= 0.0L)
        {
            outFile << datetime << "," << swh << "," << mwd << "," << pp1d;
            for (int i = 0; i < 11; i++)
                outFile << ",NaN";
            outFile << "\n";
            continue;
        }
        long double relativeDir = fmodl((mwd - coast_dir) + 360.0L, 360.0L);
        if (relativeDir < 180.0L)
        {
            outFile << datetime << "," << swh << "," << mwd << "," << pp1d;
            for (int i = 0; i < 11; i++)
                outFile << ",0.0";
            outFile << "\n";
            vector<long double> record = {swh, mwd, pp1d, 0.0L, 0.0L, 0.0L, 0.0L,
                                          0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L, 0.0L};
            for (size_t i = 0; i < NUM_COLS; i++)
                statsData[i].push_back(record[i]);
            continue;
        }
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
                << L0 << "," << depth_d << "," << L << "," << kh << ","
                << alpha_offshore_deg << "," << alpha_local_deg << ","
                << swh_local << "," << mwd_local << ","
                << Ks << "," << Kr << "," << Hb << "\n";
        vector<long double> record = {swh, mwd, pp1d, L0, depth_d, L, kh,
                                      alpha_offshore_deg, alpha_local_deg, swh_local,
                                      mwd_local, Ks, Kr, Hb};
        for (size_t i = 0; i < NUM_COLS; i++)
            statsData[i].push_back(record[i]);
    }
    inFile.close();
    outFile.close();

    // ----------------------------------------------------------------
    // Create report.txt with descriptive statistics in side-by-side format
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
    reportFile << "Descriptive Statistics Report\n";
    reportFile << string(LINE_WIDTH, '=') << "\n\n";
    for (size_t i = 0; i < NUM_COLS; i += VARIABLES_PER_ROW)
    {
        // Print variable names for this group
        for (size_t j = 0; j < VARIABLES_PER_ROW && (i + j) < NUM_COLS; ++j)
            reportFile << left << setw(colWidth) << varNames[i + j];
        reportFile << "\n";
        vector<DescriptiveStats> groupStats;
        for (size_t j = 0; j < VARIABLES_PER_ROW && (i + j) < NUM_COLS; ++j)
            groupStats.push_back(computeStats(statsData[i + j]));
        vector<string> statLabels = {"Count", "Mean", "Std. Dev.", "Min", "Q1", "Median", "Q3", "Max"};
        for (size_t statIdx = 0; statIdx < statLabels.size(); ++statIdx)
        {
            for (size_t j = 0; j < groupStats.size(); ++j)
            {
                ostringstream oss;
                if (statLabels[statIdx] == "Count")
                    oss << statLabels[statIdx] << ": " << groupStats[j].count;
                else
                {
                    oss << statLabels[statIdx] << ": " << fixed << setprecision(10);
                    if (statLabels[statIdx] == "Mean")
                        oss << groupStats[j].mean;
                    else if (statLabels[statIdx] == "Std. Dev.")
                        oss << groupStats[j].stddev;
                    else if (statLabels[statIdx] == "Min")
                        oss << groupStats[j].min;
                    else if (statLabels[statIdx] == "Q1")
                        oss << groupStats[j].Q1;
                    else if (statLabels[statIdx] == "Median")
                        oss << groupStats[j].median;
                    else if (statLabels[statIdx] == "Q3")
                        oss << groupStats[j].Q3;
                    else if (statLabels[statIdx] == "Max")
                        oss << groupStats[j].max;
                }
                reportFile << left << setw(colWidth) << oss.str();
            }
            reportFile << "\n";
        }
        reportFile << string(LINE_WIDTH, '-') << "\n";
    }
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