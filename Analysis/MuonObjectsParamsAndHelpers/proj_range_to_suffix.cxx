#include <string>
#include <sstream>
#include <iomanip>
#include <utility>
#include <stdexcept>

// function to map a projection range (pair of float) into a suffix string
std::string pairToSuffix(const std::pair<float, float>& p) {
    auto formatNum = [](float x) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2) << x;
        std::string s = oss.str();
        for (auto& c : s) if (c == '.') c = '_'; // replace '.' with '_'
        
        // replace '-' with "minus"
        size_t dash_pos;
        while ((dash_pos = s.find('-')) != std::string::npos)
            s.replace(dash_pos, 1, "minus");

        return s;
    };
    return formatNum(p.first) + "_TO_" + formatNum(p.second);
}

// function to map a suffix string (fixed format) back into a projection range (pair of float)
std::pair<float, float> suffixToPair(std::string s) {
    size_t pos = s.find("_TO_");
    if (pos == std::string::npos)
        throw std::runtime_error("suffixToPair: missing _TO_ in string: " + s);

    std::string s1 = s.substr(0, pos);
    std::string s2 = s.substr(pos + 4);

    // replace "minus" with "-" first
    size_t mpos;
    while ((mpos = s1.find("minus")) != std::string::npos)
        s1.replace(mpos, 5, "-");
    while ((mpos = s2.find("minus")) != std::string::npos)
        s2.replace(mpos, 5, "-");

    // then replace '_' with '.'
    for (auto& c : s1) if (c == '_') c = '.';
    for (auto& c : s2) if (c == '_') c = '.';

    return {std::stof(s1), std::stof(s2)};
}

// function to map a projected-histogram name (with fixed-format projection suffix) back into a projection range (pair of float)
std::pair<float, float> hProjNameToPair(std::string s) {
    // find the "_px" or "_py" marker
    size_t start = s.find("_py_");
    if (start == std::string::npos)
        start = s.find("_px_");
    if (start == std::string::npos)
        throw std::runtime_error("hProjNameToPair: no _px or _py found in string: " + s);

    // extract substring after "_px" or "_py"
    s = s.substr(start + 3);  // skip past "_px" or "_py"

    return suffixToPair(s);
}


void test_proj_range_to_suffix(){
    std::vector<std::pair<float, float>> proj_ranges = {
        {-2.45f, -2.0f}, 
        {-2.0f, -1.6f}, 
        {-1.6f, -1.3f}, 
        {-0.9f, -0.5f}, 
        {-0.5f, -0.1f}, 
        {0.1f, 0.5f}, 
        {0.5f, 1.0f}, 
        {1.3f, 1.6f}, 
        {1.6f, 2.0f}, 
        {2.0f, 225.42f}
    };

    for (auto range : proj_ranges){
        std::string suffix = pairToSuffix(range);
        std::cout << "Suffix: " << suffix << std::endl;
        std::string histProjName = "h2D_py" + suffix;
        auto rrange = hProjNameToPair(histProjName);
        std::cout << "Reverting suffix back to: " << rrange.first << ", " << rrange.second << std::endl << std::endl;

    }
    
    // suffix == "-2_4_TO_-2_0"

}
