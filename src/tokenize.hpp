#ifndef TOKENIZE_HPP_INCLUDED
#define TOKENIZE_HPP_INCLUDED

#include <string>

// grazie a https://stackoverflow.com/a/1493195/238609

template < class ContainerT >
void tokenize(const std::string& str, ContainerT& tokens,
              const std::string& delimiters = " ",
              bool trimEmpty = false) {

    std::string::size_type pos, lastPos = 0, length = str.length();

    using value_type = typename ContainerT::value_type;
    using size_type  = typename ContainerT::size_type;

    while(lastPos < length + 1) {
        pos = str.find_first_of(delimiters, lastPos);
        if(pos == std::string::npos) {
            pos = length;
        }
        if(pos != lastPos || !trimEmpty) {
            tokens.push_back(value_type(str.data()+lastPos,
                                        (size_type)pos-lastPos ));
        }
        lastPos = pos + 1;
    }
}

#endif
