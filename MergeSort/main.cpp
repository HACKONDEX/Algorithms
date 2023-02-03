#include <iostream>
#include <vector>

void merge(std::vector<int>& arr, std::vector<int>& mergedArray, int begin,
           int end, int mid) {
    int i = begin;
    int j = mid + 1;
    int t = begin;
    while (i <= mid && j <= end) {
        if (arr[i] <= arr[j])
            mergedArray[t++] = arr[i++];
        else
            mergedArray[t++] = arr[j++];
    }
    for (; i <= mid; ++i) mergedArray[t++] = arr[i];
    for (; j <= end; ++j) mergedArray[t++] = arr[j];
    for (t = begin; t <= end; ++t) arr[t] = mergedArray[t];
}

void mergeSort(std::vector<int>& arr, std::vector<int>& tmp, int begin,
               int end) {
    if (begin >= end) return;
    int mid = (begin + end) / 2;
    mergeSort(arr, tmp, begin, mid);
    mergeSort(arr, tmp, mid + 1, end);
    merge(arr, tmp, begin, end, mid);
}
int main() {
    std::vector<int> arr;
    int a = 0;
    while (std::cin >> a) {
        arr.push_back(a);
    }
    auto tmp = arr;
    mergeSort(arr, tmp, 0, arr.size() - 1);
    for (int i = 0; i < arr.size(); ++i) std::cout << arr[i] << " ";
    std::cout << std::endl;
    return 0;
}
