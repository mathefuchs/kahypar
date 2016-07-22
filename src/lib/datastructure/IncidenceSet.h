/***************************************************************************
 *  Copyright (C) 2016 Sebastian Schlag <sebastian.schlag@kit.edu>
 **************************************************************************/

#ifndef SRC_LIB_DATASTRUCTURE_INCIDENCESET_H_
#define SRC_LIB_DATASTRUCTURE_INCIDENCESET_H_

#include <x86intrin.h>

#include <algorithm>
#include <limits>
#include <utility>

#include "lib/macros.h"
#include "lib/utils/Math.h"

namespace datastructure {
template <typename T, size_t InitialSizeFactor = 2,
          T empty = std::numeric_limits<T>::max(),
          T deleted = std::numeric_limits<T>::max() - 1>
class IncidenceSet {
 private:
  using Element = std::pair<T, size_t>;
  using Position = size_t;

 public:
  explicit IncidenceSet(const T max_size) :
    _memory(nullptr),
    _size(0),
    _max_size(utils::nextPowerOfTwoCeiled(max_size + 1)),
    _max_sparse_size(InitialSizeFactor * _max_size) {
    static_assert(std::is_pod<T>::value, "T is not a POD");
    _memory = static_cast<T*>(malloc(sizeOfDense() + sizeOfSparse()));

    for (size_t i = 0; i < _max_size; ++i) {
      new(dense() + i)T(std::numeric_limits<T>::max());
    }
    // sentinel for peek() operation of incidence set
    new(dense() + _max_size)T(std::numeric_limits<T>::max());

    for (size_t i = 0; i < _max_sparse_size; ++i) {
      new (sparse() + i)Element(empty, empty);
    }
  }

  ~IncidenceSet() {
    free(_memory);
  }

  IncidenceSet(const IncidenceSet& other) = delete;
  IncidenceSet& operator= (const IncidenceSet&) = delete;

  IncidenceSet& operator= (IncidenceSet&&) = delete;

  IncidenceSet(IncidenceSet&& other) :
    _memory(other._memory),
    _size(other._size),
    _max_size(other._max_size),
    _max_sparse_size(other._max_sparse_size) {
    other._size = 0;
    other._max_size = 0;
    other._max_sparse_size = 0;
    other._memory = nullptr;
  }

  void add(const T element) {
    ASSERT(!contains(element), V(element));
    if (_size == _max_size) {
      resize();
    }

    insert(element, _size);
    dense()[_size++] = element;
  }

  void insertIfNotContained(const T element) {
    if (!contains(element)) {
      add(element);
    }
  }

  void remove(const T element) {
    swapToEnd(element);
    removeAtEnd(element);
  }

  void undoRemoval(const T element) {
    insert(element, _size);
    dense()[_size++] = element;
  }

  void swapToEnd(const T v) {
    using std::swap;
    const T index_v = get(v).second;
    swap(dense()[index_v], dense()[_size - 1]);
    update(dense()[index_v], index_v);
    update(v, _size - 1);
    ASSERT(get(v).second == _size - 1, V(v));
  }

  void removeAtEnd(const T v) {
    ASSERT(dense()[_size - 1] == v, V(v));
    remove2(v);
    // ASSERT(v_index == _size - 1, V(v));
    --_size;
  }

  // reuse position of v to store u
  void reuse(const T u, const T v) {
    ASSERT(get(v).second == _size - 1, V(v));
    const T index = remove2(v);
    ASSERT(index == _size - 1, V(index));
    insert(u, index);
    dense()[index] = u;
  }

  T peek() {
    // This works, because we ensure that 'one past the current size'
    // is an element of dense(). In case dense is full, there is
    // a sentinel in place to ensure correct behavior.
    return dense()[_size];
  }

  void undoReuse(const T u, const T v) {
    const T index = remove2(u);
    insert(v, index);
    dense()[index] = v;
  }

  void swap(IncidenceSet& other) noexcept {
    using std::swap;
    swap(_memory, other._memory);
    swap(_size, other._size);
    swap(_max_size, other._max_size);
    swap(_max_sparse_size, other._max_sparse_size);
  }

  bool contains(const T key) const {
    const Position position = find(key);
    if (position == -1 || sparse()[position].first == empty) {
      return false;
    }
    return true;
  }

  T size() const {
    return _size;
  }

  const T* begin() const {
    return dense();
  }

  const T* end() const {
    return dense() + _size;
  }

  T capacity() const {
    return _max_size;
  }

  void printAll() const {
    for (auto it = begin(); it != begin() + _size; ++it) {
      LOGVAR(*it);
    }
  }

 private:
  void insert(const T key, const size_t value) {
    ASSERT(!contains(key), V(key));
    sparse()[nextFreeSlot(key)] = { key, value };
  }

  size_t remove2(const T key) {
    ASSERT(contains(key), V(key));
    const Position position = find(key);
    sparse()[position].first = deleted;
    return sparse()[position].second;
  }

  void update(const T key, const size_t value) {
    ASSERT(contains(key), V(key));
    ASSERT(sparse()[find(key)].first == key, V(key));
    sparse()[find(key)].second = value;
  }

  const Element & get(const T& key) const {
    return sparse()[find(key)];
  }

  Position find(const T key) const {
    const Position start_position = utils::crc32(key) % _max_sparse_size;
    const Position before = start_position != 0 ? start_position - 1 : _max_size - 1;
    for (Position position = start_position; position < _max_sparse_size; position = (position + 1) % _max_sparse_size) {
      if (sparse()[position].first == empty || sparse()[position].first == key) {
        return position;
      } else if (position == before) {
        return -1;
      }
    }
    ASSERT(true == false, "This should never happen.");
  }

  Position nextFreeSlot(const T key) {
    const Position start_position = utils::crc32(key) % _max_sparse_size;
    for (Position position = start_position; position < _max_sparse_size; position = (position + 1) % _max_sparse_size) {
      if (sparse()[position].first == empty || sparse()[position].first == deleted) {
        return position;
      }
    }
    ASSERT(true == false, "This should never happen.");
  }


  size_t sizeOfDense() const {
    return (_max_size + 1  /*sentinel for peek */) * sizeof(T);
  }

  size_t sizeOfSparse() const {
    return _max_sparse_size * sizeof(std::pair<T, size_t>);
  }


  void resize() {
    IncidenceSet new_set(_max_size + 1);
    for (const auto e : * this) {
      new_set.add(e);
    }
    swap(new_set);
  }

  T* dense() {
    return _memory;
  }

  const T* dense() const {
    return _memory;
  }

  Element* sparse() {
    return reinterpret_cast<Element*>(dense() + _max_size + 1);
  }

  const Element* sparse() const {
    return reinterpret_cast<const Element*>(dense() + _max_size + 1);
  }

  T* _memory;
  T _size;
  T _max_size;
  T _max_sparse_size;
};
}  // namespace datastructure

#endif  // SRC_LIB_DATASTRUCTURE_INCIDENCESET_H_
