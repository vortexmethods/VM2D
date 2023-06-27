#ifndef LBVH_QUERY_CUH
#define LBVH_QUERY_CUH
#include "predicator.cuh"
extern int calcs;



__device__ const float machineEps = 2.0e-7f;
__device__ const float onePlusMachineEps = 1.0f + 2.0e-7f;

namespace lbvh
{
// query object indices that potentially overlaps with query aabb.
//
// requirements:
// - OutputIterator should be writable and its object_type should be uint32_t
//
template<typename Real, typename Objects, bool IsConst, typename OutputIterator>
__device__
unsigned int query_device(
        const detail::basic_device_bvh<Real, Objects, IsConst>& bvh,
        const query_overlap<Real> q, OutputIterator outiter,
        const unsigned int max_buffer_size = 0xFFFFFFFF) noexcept
{
    using bvh_type   = detail::basic_device_bvh<Real, Objects, IsConst>;
    using index_type = typename bvh_type::index_type;
    using aabb_type  = typename bvh_type::aabb_type;
    using node_type  = typename bvh_type::node_type;

    index_type  stack[64]; // is it okay?
    index_type* stack_ptr = stack;
    *stack_ptr++ = 0; // root node is always 0

    unsigned int num_found = 0;
    do
    {
        const index_type node  = *--stack_ptr;
        const index_type L_idx = bvh.nodes[node].left_idx;
        const index_type R_idx = bvh.nodes[node].right_idx;

        if(intersects(q.target, bvh.aabbs[L_idx]))
        {
            const auto obj_idx = bvh.nodes[L_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF)
            {
                if(num_found < max_buffer_size)
                {
                    *outiter++ = obj_idx;
                }
                ++num_found;
            }
            else // the node is not a leaf.
            {
                *stack_ptr++ = L_idx;
            }
        }
        if(intersects(q.target, bvh.aabbs[R_idx]))
        {
            const auto obj_idx = bvh.nodes[R_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF)
            {
                if(num_found < max_buffer_size)
                {
                    *outiter++ = obj_idx;
                }
                ++num_found;
            }
            else // the node is not a leaf.
            {
                *stack_ptr++ = R_idx;
            }
        }
    }
    while (stack < stack_ptr);
    return num_found;
}

template<typename Real, typename Objects, bool IsConst, typename OutputIterator>
__device__
unsigned int query_device2(
        const detail::basic_device_bvh2<Real, Objects, IsConst>& bvh,
        const query_overlap<Real> q, OutputIterator outiter,
        const unsigned int max_buffer_size = 0xFFFFFFFF) noexcept
{
    using bvh_type   = detail::basic_device_bvh2<Real, Objects, IsConst>;
    using index_type = typename bvh_type::index_type;
    using aabb_type  = typename bvh_type::aabb2_type;
    using node_type  = typename bvh_type::node_type;

    index_type  stack[64]; // is it okay?
    index_type* stack_ptr = stack;
    *stack_ptr++ = 0; // root node is always 0

    unsigned int num_found = 0;
    do
    {
        const index_type node  = *--stack_ptr;
        const index_type L_idx = bvh.nodes[node].left_idx;
        const index_type R_idx = bvh.nodes[node].right_idx;

        if(intersects2(q.target, bvh.aabbs2[L_idx]))
        {
            const auto obj_idx = bvh.nodes[L_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF)
            {
                if(num_found < max_buffer_size)
                {
                    *outiter++ = obj_idx;
                }
                ++num_found;
            }
            else // the node is not a leaf.
            {
                *stack_ptr++ = L_idx;
            }
        }
        if(intersects(q.target, bvh.aabbs[R_idx]))
        {
            const auto obj_idx = bvh.nodes[R_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF)
            {
                if(num_found < max_buffer_size)
                {
                    *outiter++ = obj_idx;
                }
                ++num_found;
            }
            else // the node is not a leaf.
            {
                *stack_ptr++ = R_idx;
            }
        }
    }
    while (stack < stack_ptr);
    return num_found;
}

// query object index that is the nearst to the query point.
//
// requirements:
// - DistanceCalculator must be able to calc distance between a point to an object.
//
template<typename Real, typename Objects, bool IsConst,
         typename DistanceCalculator>
__device__
thrust::pair<unsigned int, Real> query_device(
        const detail::basic_device_bvh<Real, Objects, IsConst>& bvh,
        const query_nearest<Real>& q, DistanceCalculator calc_dist) noexcept
{
    using bvh_type   = detail::basic_device_bvh<Real, Objects, IsConst>;
    using real_type  = typename bvh_type::real_type;
    using index_type = typename bvh_type::index_type;
    using aabb_type  = typename bvh_type::aabb_type;
    using node_type  = typename bvh_type::node_type;

    // pair of {node_idx, mindist}
    thrust::pair<index_type, real_type>  stack[64];
    thrust::pair<index_type, real_type>* stack_ptr = stack;
    *stack_ptr++ = thrust::make_pair(0, mindist(bvh.aabbs[0], q.target));

    unsigned int nearest = 0xFFFFFFFF;
    real_type dist_to_nearest_object = infinity<real_type>();
    do
    {
        const auto node  = *--stack_ptr;
        if(node.second > dist_to_nearest_object)
        {
            // if aabb mindist > already_found_mindist, it cannot have a nearest
            continue;
        }

        const index_type L_idx = bvh.nodes[node.first].left_idx;
        const index_type R_idx = bvh.nodes[node.first].right_idx;

        const aabb_type& L_box = bvh.aabbs[L_idx];
        const aabb_type& R_box = bvh.aabbs[R_idx];

        const real_type L_mindist = mindist(L_box, q.target);
        const real_type R_mindist = mindist(R_box, q.target);

        const real_type L_minmaxdist = minmaxdist(L_box, q.target);
        const real_type R_minmaxdist = minmaxdist(R_box, q.target);

        // there should be an object that locates within minmaxdist.

        if(L_mindist <= R_minmaxdist) // L is worth considering
        {
            const auto obj_idx = bvh.nodes[L_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF) // leaf node
            {
                const real_type dist = calc_dist(q.target, bvh.objects[obj_idx]);
                if(dist <= dist_to_nearest_object)
                {
                    dist_to_nearest_object = dist;
                    nearest = obj_idx;
                }
            }
            else
            {
                *stack_ptr++ = thrust::make_pair(L_idx, L_mindist);
            }
        }
        if(R_mindist <= L_minmaxdist) // R is worth considering
        {
            const auto obj_idx = bvh.nodes[R_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF) // leaf node
            {
                const real_type dist = calc_dist(q.target, bvh.objects[obj_idx]);
                if(dist <= dist_to_nearest_object)
                {
                    dist_to_nearest_object = dist;
                    nearest = obj_idx;
                }
            }
            else
            {
                *stack_ptr++ = thrust::make_pair(R_idx, R_mindist);
            }
        }
        assert(stack_ptr < stack + 64);
    }
    while (stack < stack_ptr);
    return thrust::make_pair(nearest, dist_to_nearest_object);
}


//__device__
//bool isCP(float2 pt)
//{
//    return ((::fabs(0.41141298 - pt.x) < 1e-6) && (::fabs(-0.01829865 - pt.y) < 1e-6));
//}

template<typename Real, typename Objects, bool IsConst,
         typename DistanceCalculator2>
__device__
thrust::pair<unsigned int, Real> query_device2(
    const detail::basic_device_bvh2<Real, Objects, IsConst>& bvh,
    const query_nearest2<Real>& q, DistanceCalculator2 calc_dist) noexcept
{
    using bvh_type = detail::basic_device_bvh2<Real, Objects, IsConst>;
    using real_type = typename bvh_type::real_type;
    using index_type = typename bvh_type::index_type;
    using aabb_type = typename bvh_type::aabb2_type;
    using node_type = typename bvh_type::node_type;

    // pair of {node_idx, mindist}
    thrust::pair<index_type, real_type>  stack[64];
    thrust::pair<index_type, real_type>* stack_ptr = stack;
    *stack_ptr++ = thrust::make_pair(0, mindist(bvh.aabbs2[0], q.target));

    //if (isCP(q.target))
    //    printf("FOUND!\n");

    unsigned int nearest = 0xFFFFFFFF;
    real_type dist_to_nearest_object = infinity<real_type>();
    int location;

    //int counter = 0;

    do
    {
        //if (isCP(q.target))
        //{
        //    ++counter;
        //    printf("counter = %d\n", counter);
        //}

        const auto node = *--stack_ptr;
        if (node.second > dist_to_nearest_object)
        {
            // if aabb mindist > already_found_mindist, it cannot have a nearest
            continue;
        }

        const index_type L_idx = bvh.nodes[node.first].left_idx;
        const index_type R_idx = bvh.nodes[node.first].right_idx;

        const aabb_type& L_box = bvh.aabbs2[L_idx];
        const aabb_type& R_box = bvh.aabbs2[R_idx];

        //if (isCP(q.target))
        //{
        //    printf(" leftChild = { {%f, %f}, {%f, %f} }\n", L_box.lower.x, L_box.lower.y, L_box.upper.x, L_box.upper.y);
        //    printf("rightChild = { {%f, %f}, {%f, %f} }\n", R_box.lower.x, R_box.lower.y, R_box.upper.x, R_box.upper.y);
        //}

        const real_type L_mindist = mindist(L_box, q.target);
        const real_type R_mindist = mindist(R_box, q.target);

        const real_type L_minmaxdist = minmaxdist2(L_box, q.target);
        const real_type R_minmaxdist = minmaxdist2(R_box, q.target);

        // there should be an object that locates within minmaxdist.

        //if (isCP(q.target))
        //{
        //    printf("L_mindist = %.16f, R_minmaxdist = %.16f\n", L_mindist, R_minmaxdist);
        //}
        //if(L_mindist <= R_minmaxdist) // L is worth considering
        //if ((L_mindist < R_minmaxdist) || (::fabsf(L_mindist - R_minmaxdist) <= machineEps * ::fmaxf(L_mindist, R_minmaxdist)))   //  << todo
        if (L_mindist <= R_minmaxdist * onePlusMachineEps) // L is worth considering
        {

            const auto obj_idx = bvh.nodes[L_idx].object_idx;
            if (obj_idx != 0xFFFFFFFF) // leaf node
            {
                //if (isCP(q.target))
                //{
                //    printf("Left -> panel\n");
                //}
                const thrust::pair<real_type, int> dist = calc_dist(q.target, bvh.objects[obj_idx]);
                if (dist.first <= dist_to_nearest_object)
                {
                    dist_to_nearest_object = dist.first;
                    location = dist.second;
                    nearest = obj_idx;
                }
            }
            else
            {
                //if (isCP(q.target))
                //{
                //    printf("Left -> continue\n");
                //}
                *stack_ptr++ = thrust::make_pair(L_idx, L_mindist);
            }
        }


        //if (isCP(q.target))
        //{
        //    printf("R_mindist = %.16f, L_minmaxdist = %.16f\n", R_mindist, L_minmaxdist);
        //}

        //if(R_mindist <= L_minmaxdist) // R is worth considering
        //if ((R_mindist < L_minmaxdist) || (::fabsf(R_mindist - L_minmaxdist) <= machineEps * ::fmaxf(R_mindist, L_minmaxdist)))      // << todo
        if(R_mindist <= L_minmaxdist * onePlusMachineEps) // R is worth considering
        {
            const auto obj_idx = bvh.nodes[R_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF) // leaf node
            {
                //if (isCP(q.target))
                //{
                //    printf("Right -> panel\n");
                //}
                
                const thrust::pair<real_type, int> dist = calc_dist(q.target, bvh.objects[obj_idx]);
                if(dist.first <= dist_to_nearest_object)
                {
                    dist_to_nearest_object = dist.first;
                    location = dist.second;
                    nearest = obj_idx;
                }
            }
            else
            {
                //if (isCP(q.target))
                //{
                //    printf("Right -> continue\n");
                //}
                *stack_ptr++ = thrust::make_pair(R_idx, R_mindist);
            }
        }
        assert(stack_ptr < stack + 64);
    }
    while (stack < stack_ptr);
    return thrust::make_pair(nearest, location);
}




template<typename Real, typename Objects, bool IsConst,
    typename HitChecker2, typename HitChecker2node>
__device__
unsigned int query_device2ray(
    const detail::basic_device_bvh2<Real, Objects, IsConst>& bvh,
    const query_nearest2ray<Real>& q, HitChecker2 check_hit_obj, HitChecker2node check_hit_node) noexcept
{
    using bvh_type = detail::basic_device_bvh2<Real, Objects, IsConst>;
    using real_type = typename bvh_type::real_type;
    using index_type = typename bvh_type::index_type;
    using aabb_type = typename bvh_type::aabb2_type;
    using node_type = typename bvh_type::node_type;

    // pair of {node_idx, mindist}
    thrust::pair<index_type, real_type>  stack[256];
    thrust::pair<index_type, real_type>* stack_ptr = stack;

    unsigned int hitted = 0xFFFFFFFF;
    real_type dist_to_hitted_object = infinity<real_type>();

    //auto distToRoot = check_hit_node(q.start, q.finish, bvh.aabbs2[0]);
    //if (distToRoot < infinity<real_type>())
    //    *stack_ptr++ = thrust::make_pair(0, distToRoot);

    bool crossRoot = aabbContainsSegment(q.start, q.finish, bvh.aabbs2[0]);
    if (crossRoot)
        *stack_ptr++ = thrust::make_pair(0, 0.0);
    
    while (stack < stack_ptr)
    {
        const auto node = *--stack_ptr;        

        //printf("node.second = %f, dist = %f\n", node.second, dist_to_hitted_object);
        if (node.second > dist_to_hitted_object)
        {
            // if aabb mindist > already_found_mindist, it cannot have a nearest
            continue;
        }             
          
        const auto obj_idx = bvh.nodes[node.first].object_idx;
        //printf("C, node.first = %d, obj[node.first] = %d\n", (int)node.first, (int)obj_idx);

        if (obj_idx != 0xFFFFFFFF)
        {
            //printf("Leaf_processing\n");
            real_type hit = check_hit_obj(q.start, q.finish, bvh.objects[obj_idx]);      ///////// check_hit_obj - проверка пересечения лучом отрезка
            //printf("hit = %f\n", hit);
            if (hit < dist_to_hitted_object)
            {                
                dist_to_hitted_object = hit;                
                hitted = obj_idx;
            }
        }
        else
        {
            //printf("F\n");
            const index_type L_idx = bvh.nodes[node.first].left_idx;
            const index_type R_idx = bvh.nodes[node.first].right_idx;

            const aabb_type& L_box = bvh.aabbs2[L_idx];
            const aabb_type& R_box = bvh.aabbs2[R_idx];

            real_type h_L = check_hit_node(q.start, q.finish, L_box);  ///////// check_hit_node - проверка пересечения лучом прямоугольника
            real_type h_R = check_hit_node(q.start, q.finish, R_box); 

            //printf("q_start = {%f, %f}, q_finish = {%f, %f}\n", q.start.x, q.start.y, q.finish.x, q.finish.y);
            //printf("L_box = { {%f, %f}, {%f, %f} }\n", L_box.lower.x, L_box.lower.y, L_box.upper.x, L_box.upper.y);
            //printf("R_box = { {%f, %f}, {%f, %f} }\n", R_box.lower.x, R_box.lower.y, R_box.upper.x, R_box.upper.y);
            //printf("h_L = %f, h_R = %f\n", (float)h_L, (float)h_R);
            
            if ((h_L) == infinity<real_type>() && (h_R) == infinity<real_type>())
            {
            }
            else if ((h_L) < infinity<real_type>() && (h_R) == infinity<real_type>())
            { 
                *stack_ptr++ = thrust::make_pair(L_idx, h_L);
                //printf("Left_added\n");
            }
            else if ((h_R) < infinity<real_type>() && (h_L) == infinity<real_type>())
            { 
                *stack_ptr++ = thrust::make_pair(R_idx, h_R);
                //printf("Right_added\n");
            }
            else
            {
                if (h_L > h_R)
                {                    
                    *stack_ptr++ = thrust::make_pair(L_idx, h_L);
                    *stack_ptr++ = thrust::make_pair(R_idx, h_R);
                    //printf("Left_added\n");
                    //printf("Right_added\n");
                }
                else
                {                    
                    *stack_ptr++ = thrust::make_pair(R_idx, h_R);
                    *stack_ptr++ = thrust::make_pair(L_idx, h_L);
                    //printf("Right_added\n");
                    //printf("Left_added\n");                    
                }
            }

        }

        assert(stack_ptr < stack + 256);
    };
    //printf("hitted = %d\n", (int)hitted);
    return hitted;
}





template<typename Real, typename Objects, typename AABBGetter,
         typename MortonCodeCalculator, typename OutputIterator>
unsigned int query_host(
    const bvh<Real, Objects, AABBGetter, MortonCodeCalculator>& tree,
    const query_overlap<Real> q, OutputIterator outiter,
    const unsigned int max_buffer_size = 0xFFFFFFFF)
{
    using bvh_type   = ::lbvh::bvh<Real, Objects, AABBGetter, MortonCodeCalculator>;
    using index_type = typename bvh_type::index_type;
    using aabb_type  = typename bvh_type::aabb_type;
    using node_type  = typename bvh_type::node_type;

    if(!tree.query_host_enabled())
    {
        throw std::runtime_error("lbvh::bvh query_host is not enabled");
    }

    std::vector<std::size_t> stack;
    stack.reserve(64);
    stack.push_back(0);

    unsigned int num_found = 0;
    do
    {
        const index_type node  = stack.back(); stack.pop_back();
        const index_type L_idx = tree.nodes_host()[node].left_idx;
        const index_type R_idx = tree.nodes_host()[node].right_idx;

        if(intersects(q.target, tree.aabbs_host()[L_idx]))
        {
            const auto obj_idx = tree.nodes_host()[L_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF)
            {
                if(num_found < max_buffer_size)
                {
                    *outiter++ = obj_idx;
                }
                ++num_found;
            }
            else // the node is not a leaf.
            {
                stack.push_back(L_idx);
            }
        }
        if(intersects(q.target, tree.aabbs_host()[R_idx]))
        {
            const auto obj_idx = tree.nodes_host()[R_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF)
            {
                if(num_found < max_buffer_size)
                {
                    *outiter++ = obj_idx;
                }
                ++num_found;
            }
            else // the node is not a leaf.
            {
                stack.push_back(R_idx);
            }
        }
    }
    while (!stack.empty());
    return num_found;
}


template<typename Real, typename Objects, typename AABB2Getter,
         typename MortonCodeCalculator2, typename OutputIterator>
unsigned int query_host2(
    const bvh2<Real, Objects, AABB2Getter, MortonCodeCalculator2>& tree,
    const query_overlap2<Real> q, OutputIterator outiter,
    const unsigned int max_buffer_size = 0xFFFFFFFF)
{
    using bvh_type   = ::lbvh::bvh2<Real, Objects, AABB2Getter, MortonCodeCalculator2>;
    using index_type = typename bvh_type::index_type;
    using aabb2_type  = typename bvh_type::aabb2_type;
    using node_type  = typename bvh_type::node_type;

    if(!tree.query_host_enabled())
    {
        throw std::runtime_error("lbvh::bvh query_host is not enabled");
    }

    std::vector<std::size_t> stack;
    stack.reserve(64);
    stack.push_back(0);

    unsigned int num_found = 0;
    do
    {
        const index_type node  = stack.back(); stack.pop_back();
        const index_type L_idx = tree.nodes_host()[node].left_idx;
        const index_type R_idx = tree.nodes_host()[node].right_idx;

        if(intersects2(q.target, tree.aabbs2_host()[L_idx]))
        {
            const auto obj_idx = tree.nodes_host()[L_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF)
            {
                if(num_found < max_buffer_size)
                {
                    *outiter++ = obj_idx;
                }
                ++num_found;
            }
            else // the node is not a leaf.
            {
                stack.push_back(L_idx);
            }
        }
        if(intersects2(q.target, tree.aabbs2_host()[R_idx]))
        {
            const auto obj_idx = tree.nodes_host()[R_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF)
            {
                if(num_found < max_buffer_size)
                {
                    *outiter++ = obj_idx;
                }
                ++num_found;
            }
            else // the node is not a leaf.
            {
                stack.push_back(R_idx);
            }
        }
    }
    while (!stack.empty());
    return num_found;
}


template<typename Real, typename Objects, typename AABBGetter,
         typename MortonCodeCalculator, typename DistanceCalculator>
std::pair<unsigned int, Real> query_host(
        const bvh<Real, Objects, AABBGetter, MortonCodeCalculator>& tree,
        const query_nearest<Real>& q, DistanceCalculator calc_dist) noexcept
{
    using bvh_type   = ::lbvh::bvh<Real, Objects, AABBGetter, MortonCodeCalculator>;
    using real_type  = typename bvh_type::real_type;
    using index_type = typename bvh_type::index_type;
    using aabb_type  = typename bvh_type::aabb_type;
//    using node_type  = typename bvh_type::node_type;

    if(!tree.query_host_enabled())
    {
        throw std::runtime_error("lbvh::bvh query_host is not enabled");
    }

    // pair of {node_idx, mindist}
    std::vector<std::pair<index_type, real_type>> stack = {
        {0, mindist(tree.aabbs_host()[0], q.target)}
    };
    stack.reserve(64);

    unsigned int nearest = 0xFFFFFFFF;
    real_type current_nearest_dist = infinity<real_type>();
    do
    {
        const auto node = stack.back(); stack.pop_back();
        if(node.second > current_nearest_dist)
        {
            // if aabb mindist > already_found_mindist, it cannot have a nearest
            continue;
        }

        const index_type L_idx = tree.nodes_host()[node.first].left_idx;
        const index_type R_idx = tree.nodes_host()[node.first].right_idx;

        const aabb_type& L_box = tree.aabbs_host()[L_idx];
        const aabb_type& R_box = tree.aabbs_host()[R_idx];

        const real_type L_mindist = mindist(L_box, q.target);
        const real_type R_mindist = mindist(R_box, q.target);

        const real_type L_minmaxdist = minmaxdist(L_box, q.target);
        const real_type R_minmaxdist = minmaxdist(R_box, q.target);

       // there should be an object that locates within minmaxdist.

        if(L_mindist <= R_minmaxdist) // L is worth considering
        {
            const auto obj_idx = tree.nodes_host()[L_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF) // leaf node
            {
                const real_type dist = calc_dist(q.target, tree.objects_host()[obj_idx]);
                if(dist <= current_nearest_dist)
                {
                    current_nearest_dist = dist;
                    nearest = obj_idx;
                }
            }
            else
            {
                stack.emplace_back(L_idx, L_mindist);
            }
        }
        if(R_mindist <= L_minmaxdist) // R is worth considering
        {
            const auto obj_idx = tree.nodes_host()[R_idx].object_idx;
            if(obj_idx != 0xFFFFFFFF) // leaf node
            {
                const real_type dist = calc_dist(q.target, tree.objects_host()[obj_idx]);
                if(dist <= current_nearest_dist)
                {
                    current_nearest_dist = dist;
                    nearest = obj_idx;
                }
            }
            else
            {
                stack.emplace_back(R_idx, R_mindist);
            }
        }
    }
    while (!stack.empty());
    return std::make_pair(nearest, current_nearest_dist);
}



} // lbvh
#endif// LBVH_QUERY_CUH
