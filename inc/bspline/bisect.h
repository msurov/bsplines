
/**
 * @brief Find index i s.t. xarr[i] <= x < xarr[i+1]
 *  @param arr is a sorted array
 *  @param x argument
 */
template <typename Array, typename T>
inline int bisect(Array const& arr, T x)
{
  if (arr.size() == 0)
    return -1;

  int i1 = 0;
  int i2 = static_cast<int>(arr.size());
  while (i2 - i1 > 1) {
    int i = (i1 + i2) / 2;
    if (arr[i] <= x) {
      i1 = i;
    } else {
      i2 = i;
    }
  }

  return i1;
}
