
#ifndef _ITERATOR_H_
#define _ITERATOR_H_

template <typename CoeffIt, typename VarIt>
class Iterator {
    private:
        CoeffIt coeffIt;
        VarIt varIt;
    public:
        Iterator(CoeffIt coeffIt, VarIt varIt)
            : coeffIt(coeffIt), varIt(varIt) {}

        // dereference操作でpairを返す
        auto operator*() const {
            return std::make_pair(*coeffIt, *varIt);
        }

        // 前置インクリメント
        Iterator& operator++() {
            ++coeffIt;
            ++varIt;
            return *this;
        }

        // イテレータの比較
        bool operator!=(const Iterator& other) const {
            return coeffIt != other.coeffIt || varIt != other.varIt;
        }
};

#endif 