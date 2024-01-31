#include<set>
#include<vector>
#include<map>
#include<string>
using namespace std; 
namespace CAE
{
	
	class SurfaceNode
	{
		class face {
		public:
			face(int a, int b, int c) {
				set<int> sort = { a,b,c };
				auto v = sort.begin();
				this->a = *v++;
				this->b = *v++;
				this->c = *v;
			}
			int a; 
			int b; 
			int c;
			bool operator ==(const face& f) const {
				if (a == f.a && b == f.b && c == f.c) {
					return true;
				}
				return false;
			}
			bool operator <(const face& f) const{
				if (a == f.a) {
					if (b == f.b) {
						if (c == f.c) {
							return false;
						}
						else if (c > f.c) {
							return false;
						}
						return true;
					}
					else if (b > f.b) {
						return false;
					}
					return true;
				}
				else if (a > f.a) {
					return false;
				}
				return true;
			}
			bool operator <=(const face& f)const {
				if (f == *this) {
					return true;
				}
				if (*this < f) {
					return true;
				}
				return false;
			}
			bool operator >(const face& f)const {
				if (f == *this) {
					return false;
				}
				if (*this < f) {
					return false;
				}
				return true;
			}
			bool operator >=(const face& f)const {
				if (f == *this) {
					return true;
				}
				if (*this < f) {
					return false;
				}
				return true;
			}
		};
	public:
		SurfaceNode();
		~SurfaceNode();

		set<face> T3Surface;


		void insertFace(int a, int b, int c);

		bool isSurface(int a, int b, int c);

		set<int> node();

		vector<int> currentFaceNode();
	};
}