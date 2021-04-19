#ifndef TABLE_H__
#define TABLE_H__


struct Table {
  // get ith row of the table; max size of a table row is given by std::size_t
  virtual std::size_t operator()(std::size_t i) const = 0;
  virtual ~Table() { };
};


#endif
