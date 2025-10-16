#include <string>
#include <sstream>
#include <iomanip>
#include <utility>

std::string pairToSuffix(const std::pair<float, float>& p) {
    auto formatNum = [](float x) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2) << x;
        std::string s = oss.str();
        for (auto& c : s) if (c == '.') c = '_'; // replace '.' with '_'
        return s;
    };
    return formatNum(p.first) + "_TO_" + formatNum(p.second);
}

std::pair<float, float> suffixToPair(std::string s) {
    // Replace '_' back to '.' except in "_TO_"
    size_t pos = s.find("_TO_");
    std::string s1 = s.substr(0, pos);
    std::string s2 = s.substr(pos + 4);

    for (auto& c : s1) if (c == '_') c = '.';
    for (auto& c : s2) if (c == '_') c = '.';

    return {std::stof(s1), std::stof(s2)};
}

void test(){
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
        auto rrange = suffixToPair(suffix);
        std::cout << "Reverting suffix back to: " << rrange.first << ", " << rrange.second << std::endl << std::endl;

    }
    
    // suffix == "-2_4_TO_-2_0"

}
