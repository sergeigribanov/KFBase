#ifndef _KFBASE_OPT_NAMED_HPP_
#define _KFBASE_OPT_NAMED_HPP_
#include <string>

namespace kfbase {
  namespace newtonian_opt {
    class Named {
    public:
      explicit Named(const std::string&);
      virtual ~Named();
      const std::string& getName() const;
    private:
      std::string name_;
    };
  } // namespace newtonian_opt
} // namespace kfbase

#endif
