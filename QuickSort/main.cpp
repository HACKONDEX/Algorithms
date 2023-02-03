#include <iostream>
#include <vector>

void InsertSort(std::vector<int>& arr, int begin, int end) {
    int tmp = 0;
    for (int i = begin; i <= end; ++i) {
        tmp = arr[i];
        int j = i - 1;
        for (; j >= begin; --j) {
            if (tmp <= arr[j]) {
                arr[j + 1] = arr[j];
            } else {
                arr[j + 1] = tmp;
                break;
            }
        }
        if (j < begin) {
            arr[j + 1] = tmp;
        }
    }
}

int ChoosePivotMedian(int begin, int end) { return ((begin + end) / 2); }

int ChoosePivotMedianOfThree(std::vector<int>& arr, int begin, int end) {
    int index = begin;
    int mid = (end + begin) / 2;
    if (arr[begin] <= arr[mid]) {
        if (arr[mid] <= arr[end])
            return mid;
        else {
            if (arr[begin] <= arr[end])
                return end;
            else
                return begin;
        }
    } else {
        if (arr[begin] <= arr[end])
            return begin;
        else {
            if (arr[mid] <= arr[end])
                return end;
            else
                return mid;
        }
    }
}

int ChoosePivotRandom(int begin, int end) {
    return (rand() % (end - begin + 1) + begin);
}

int PartitionFromBegin(std::vector<int>& arr, int begin, int end) {
    int pivotIndex = ChoosePivotRandom(begin, end);
    int pivot = arr[pivotIndex];
    std::swap(arr[pivotIndex], arr[end]);
    int i = begin;
    int j = begin;

    while (j != end) {
        if (arr[j] <= pivot) {
            std::swap(arr[i], arr[j]);
            ++i;
        }
        ++j;
    }

    return (i - 1);
}

int PartitionFromBothEnds(std::vector<int>& arr, int begin, int end) {
    int pivotIndex = ChoosePivotMedianOfThree(arr, begin, end);
    int pivot = arr[pivotIndex];
    int i = begin, j = end;
    while (i < j) {
        while (arr[i] < pivot) ++i;
        while (arr[j] > pivot) --j;
        if (i < j) std::swap(arr[i++], arr[j--]);
    }
    return j;
}

void QuickSort(std::vector<int>& arr, int begin, int end) {
    if (begin >= end) return;
    if (end - begin + 1 <= 25) {
        InsertSort(arr, begin, end);
        return;
    }
    int partitionIndex;
    partitionIndex = PartitionFromBothEnds(arr, begin, end);
    QuickSort(arr, begin, partitionIndex);
    QuickSort(arr, partitionIndex + 1, end);
}
