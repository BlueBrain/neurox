#pragma once

class IDataLoader {
public:
  virtual void loadData() = 0; ///>Loads data
  void clearData(); ///>Clears data
};
