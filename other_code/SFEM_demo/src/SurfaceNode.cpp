#include "include/SurfaceNode.h"

namespace CAE
{
	SurfaceNode::SurfaceNode()
	{

	}

	SurfaceNode::~SurfaceNode()
	{

	}

	void SurfaceNode::insertFace(int a, int b,int c)
	{
		face f = face(a, b, c);
		auto target = T3Surface.find(f);
		if (target == T3Surface.end()) {
			T3Surface.insert(f);
		}
		else {
			T3Surface.erase(target);
		}	
	}

	bool SurfaceNode::isSurface(int a, int b, int c)
	{
		face f = face(a, b, c);
		auto target = T3Surface.find(f);
		if (target == T3Surface.end()) {
			return false;
		}
		else {
			return true;
		}
	}

	set<int> SurfaceNode::node()
	{
		set<int> node;
		for (auto i : T3Surface) {
			node.insert(i.a);
			node.insert(i.b);
			node.insert(i.c);
		}
		return node;
	}

	vector<int> SurfaceNode::currentFaceNode()
	{
		return vector<int>();
	}
}