#pragma once

// ---------- try catch for at() ----------
template <class Map>
auto& map_at_checked(Map& m, const std::string& key, const char* where)
{
    try {
        return m.at(key);
    } catch (const std::out_of_range& e) {
        throw std::out_of_range(
            std::string("map_at_checked: missing key '") + key +
            "' in " + where + " (" + e.what() + ")"
        );
    }
}

// could apply to Map/key or vector/index pairs
template <class Container, class Index>
decltype(auto) at_checked(const Container& c, const Index& idx, const char* where)
{
    try {
        return c.at(idx);
    } catch (const std::out_of_range& e) {
        throw std::out_of_range(
            std::string("at_checked: out_of_range in ") +
            where + " (" + e.what() + ")"
        );
    }
}

template <typename Map, typename Key, typename Vec>
void InsertOrAppend(Map& m, const Key& key, const Vec& v_new)
{
    auto [it, inserted] = m.try_emplace(key, v_new);

    if (!inserted) {
        it->second.insert(it->second.end(),
                          v_new.begin(),
                          v_new.end());
    }
}
