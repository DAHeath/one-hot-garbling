#ifndef TABLE_H__
#define TABLE_H__


struct Table {
  virtual bool operator()(std::size_t i, std::size_t j) const = 0;
  virtual ~Table() { };
};


#endif
