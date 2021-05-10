#ifndef _KDTREE_H_
#define _KDTREE_H_

#ifdef __cplusplus
extern "C" {
#endif

struct kdtree;
struct kdres;


/*创建kd树 */
struct kdtree *kd_create(int k);

/* 释放kd树 */
void kd_free(struct kdtree *tree);

/* r移除树中所有元素 */
void kd_clear(struct kdtree *tree);

/* 如果使用第二个非空参数调用，则从树中删除节点时，将在数据指针上调用提供的函数（请参见kd_insert）*/
void kd_data_destructor(struct kdtree *tree, void (*destr)(void*));

/* 插入节点，指定其位置和可选数据 */
int kd_insert(struct kdtree *tree, const double *pos, void *data);
int kd_insertf(struct kdtree *tree, const float *pos, void *data);
int kd_insert3(struct kdtree *tree, double x, double y, double z, void *data);
int kd_insert3f(struct kdtree *tree, float x, float y, float z, void *data);

/* 从给定点找到最近的节点。
 * 此函数返回一个最多包含一个元素的结果集的指针。*/
struct kdres *kd_nearest(struct kdtree *tree, const double *pos);
struct kdres *kd_nearestf(struct kdtree *tree, const float *pos);
struct kdres *kd_nearest3(struct kdtree *tree, double x, double y, double z);
struct kdres *kd_nearest3f(struct kdtree *tree, float x, float y, float z);

/* 从给定点找到N个最近的节点。
 * 此函数返回一个最多包含N个元素的结果集的指针，可以通过kd_res_ *函数进行操作。
 * 返回的指针可以为null，以指示错误。 除此以外，始终返回可能包含0个或多个元素的有效结果集。
 * 使用后必须用kd_res_free释放结果集。
*/
/*
struct kdres *kd_nearest_n(struct kdtree *tree, const double *pos, int num);
struct kdres *kd_nearest_nf(struct kdtree *tree, const float *pos, int num);
struct kdres *kd_nearest_n3(struct kdtree *tree, double x, double y, double z);
struct kdres *kd_nearest_n3f(struct kdtree *tree, float x, float y, float z);
*/

/* 查找范围内给定点的任何最近节点。
 * 此函数返回一个指向结果集的指针，可以通过kd_res_ *函数对它进行操作。
 * 返回的指针可以为null，以指示错误。 除此以外，始终返回可能包含0个或多个元素的有效结果集。
 * 使用后必须用kd_res_free释放结果集。
 */
struct kdres *kd_nearest_range(struct kdtree *tree, const double *pos, double range);
struct kdres *kd_nearest_rangef(struct kdtree *tree, const float *pos, float range);
struct kdres *kd_nearest_range3(struct kdtree *tree, double x, double y, double z, double range);
struct kdres *kd_nearest_range3f(struct kdtree *tree, float x, float y, float z, float range);

/* 释放由kd_nearest_range（）返回的结果集 */
void kd_res_free(struct kdres *set);

/* 返回结果集的大小（以元素为单位） */
int kd_res_size(struct kdres *set);

/* 倒退结果集迭代器 */
void kd_res_rewind(struct kdres *set);

/* 如果集合迭代器在最后一个元素之后到达结尾，则返回非零 */
int kd_res_end(struct kdres *set);

/* 推进结果集迭代器，成功返回非零，如果结果集中没有元素，则返回零。*/
int kd_res_next(struct kdres *set);

/* 返回当前结果集项的数据指针（可以为null），如果不为null，则可以选择将其位置设置为指针。*/
void *kd_res_item(struct kdres *set, double *pos);
void *kd_res_itemf(struct kdres *set, float *pos);
void *kd_res_item3(struct kdres *set, double *x, double *y, double *z);
void *kd_res_item3f(struct kdres *set, float *x, float *y, float *z);

/* 等效于kd_res_item（set，0） */
void *kd_res_item_data(struct kdres *set);


#ifdef __cplusplus
}
#endif

#endif
