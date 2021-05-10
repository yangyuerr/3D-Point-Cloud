#ifndef _KDTREE_H_
#define _KDTREE_H_

#ifdef __cplusplus
extern "C" {
#endif

struct kdtree;
struct kdres;


/*����kd�� */
struct kdtree *kd_create(int k);

/* �ͷ�kd�� */
void kd_free(struct kdtree *tree);

/* r�Ƴ���������Ԫ�� */
void kd_clear(struct kdtree *tree);

/* ���ʹ�õڶ����ǿղ������ã��������ɾ���ڵ�ʱ����������ָ���ϵ����ṩ�ĺ�������μ�kd_insert��*/
void kd_data_destructor(struct kdtree *tree, void (*destr)(void*));

/* ����ڵ㣬ָ����λ�úͿ�ѡ���� */
int kd_insert(struct kdtree *tree, const double *pos, void *data);
int kd_insertf(struct kdtree *tree, const float *pos, void *data);
int kd_insert3(struct kdtree *tree, double x, double y, double z, void *data);
int kd_insert3f(struct kdtree *tree, float x, float y, float z, void *data);

/* �Ӹ������ҵ�����Ľڵ㡣
 * �˺�������һ��������һ��Ԫ�صĽ������ָ�롣*/
struct kdres *kd_nearest(struct kdtree *tree, const double *pos);
struct kdres *kd_nearestf(struct kdtree *tree, const float *pos);
struct kdres *kd_nearest3(struct kdtree *tree, double x, double y, double z);
struct kdres *kd_nearest3f(struct kdtree *tree, float x, float y, float z);

/* �Ӹ������ҵ�N������Ľڵ㡣
 * �˺�������һ��������N��Ԫ�صĽ������ָ�룬����ͨ��kd_res_ *�������в�����
 * ���ص�ָ�����Ϊnull����ָʾ���� �������⣬ʼ�շ��ؿ��ܰ���0������Ԫ�ص���Ч�������
 * ʹ�ú������kd_res_free�ͷŽ������
*/
/*
struct kdres *kd_nearest_n(struct kdtree *tree, const double *pos, int num);
struct kdres *kd_nearest_nf(struct kdtree *tree, const float *pos, int num);
struct kdres *kd_nearest_n3(struct kdtree *tree, double x, double y, double z);
struct kdres *kd_nearest_n3f(struct kdtree *tree, float x, float y, float z);
*/

/* ���ҷ�Χ�ڸ�������κ�����ڵ㡣
 * �˺�������һ��ָ��������ָ�룬����ͨ��kd_res_ *�����������в�����
 * ���ص�ָ�����Ϊnull����ָʾ���� �������⣬ʼ�շ��ؿ��ܰ���0������Ԫ�ص���Ч�������
 * ʹ�ú������kd_res_free�ͷŽ������
 */
struct kdres *kd_nearest_range(struct kdtree *tree, const double *pos, double range);
struct kdres *kd_nearest_rangef(struct kdtree *tree, const float *pos, float range);
struct kdres *kd_nearest_range3(struct kdtree *tree, double x, double y, double z, double range);
struct kdres *kd_nearest_range3f(struct kdtree *tree, float x, float y, float z, float range);

/* �ͷ���kd_nearest_range�������صĽ���� */
void kd_res_free(struct kdres *set);

/* ���ؽ�����Ĵ�С����Ԫ��Ϊ��λ�� */
int kd_res_size(struct kdres *set);

/* ���˽���������� */
void kd_res_rewind(struct kdres *set);

/* ������ϵ����������һ��Ԫ��֮�󵽴��β���򷵻ط��� */
int kd_res_end(struct kdres *set);

/* �ƽ���������������ɹ����ط��㣬����������û��Ԫ�أ��򷵻��㡣*/
int kd_res_next(struct kdres *set);

/* ���ص�ǰ������������ָ�루����Ϊnull���������Ϊnull�������ѡ����λ������Ϊָ�롣*/
void *kd_res_item(struct kdres *set, double *pos);
void *kd_res_itemf(struct kdres *set, float *pos);
void *kd_res_item3(struct kdres *set, double *x, double *y, double *z);
void *kd_res_item3f(struct kdres *set, float *x, float *y, float *z);

/* ��Ч��kd_res_item��set��0�� */
void *kd_res_item_data(struct kdres *set);


#ifdef __cplusplus
}
#endif

#endif
