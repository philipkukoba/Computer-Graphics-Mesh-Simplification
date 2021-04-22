#pragma once

template<
    class T,
    class Container = std::vector<T>,
    class Compare = std::less<typename Container::value_type>
>
class my_priority_queue : public std::priority_queue<T, Container, Compare>
{
public:

    void remove(const T& value) 
    {
        auto it = std::find(this->c.begin(), this->c.end(), value);
        if (it != this->c.end())  // if found
            this->c.erase(it);
    }
};