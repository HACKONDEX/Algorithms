#include <iostream>
#include <vector>

template <class T>
int32_t Partition(std::vector<T>& nums, int32_t left, int32_t right) {
    int32_t pivot = nums[right];
    int32_t i = left;

    for (int32_t j = left; j <= right - 1; ++j) {
        if (nums[j] < pivot) {
            std::swap(nums[i++], nums[j]);
        }
    }

    int32_t pivot_index = i;
    std::swap(nums[right], nums[pivot_index]);
    return pivot_index;
}

template <class T>
int32_t QuickSelectPrivate(std::vector<T>& nums, int32_t left, int32_t right,
                    int32_t k) {
    if (right == left) {
        return left;
    }

    int32_t pivot_index = Partition(nums, left, right);
    int32_t i = pivot_index;
    while (i >= left && nums[i] == nums[pivot_index]) {
        --i;
    }

    int32_t less_count = i < left ? 0 : i - left + 1;
    int32_t equal_count = (pivot_index - left + 1) - less_count;

    if (k <= less_count) {
        return QuickSelectPrivate(nums, left, left + (less_count - 1), k);
    }

    if (k <= less_count + equal_count) {
        return pivot_index;
    }

    return QuickSelectPrivate(nums, pivot_index + 1, right,
                       k - less_count - equal_count);
}

template <class T>
T QuickSelect(std::vector<T>& array, int32_t k_th) {
    return array[QuickSelectPrivate(array, 0, array.size() - 1, k_th)];
}

int main() {
    // Use example
    std::vector<int64_t> array = {1, 2, 10, 9, 5, 7, 4, 3, 8, 6};
    int64_t value = QuickSelect(array, 7);

    std::cout << value << '\n'; // expect 7 as it is 7th value in sorted array
    return 0;
}