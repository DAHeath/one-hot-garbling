#ifndef MEASURE_LINK_H__
#define MEASURE_LINK_H__


#include "link.h"

template <typename L>
struct MeasureLink : public Link {
public:
  MeasureLink();
  MeasureLink(L* under) : under(under) { }

  void send(std::span<const std::byte> s) {
    n += s.size();
    under->send(s);
  }
  void recv(std::span<std::byte> s) {
    n += s.size();
    under->recv(s);
  }

  void flush() {
    under->flush();
  }

  std::size_t count() const { return n; }

  void reset_count() { n = 0; }

private:
  std::size_t n = 0;
  L* under;
};

#endif
