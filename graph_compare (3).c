/*
 * 과제 - 07: 인접 행렬(Adjacency Matrix)과 인접 리스트(Adjacency List)의 성능 비교
 * 무방향 단순 그래프, 정점 100개 (번호 0 ~ 99)
 *   - 희소 그래프(Sparse): 간선 100개
 *   - 밀집 그래프(Dense) : 간선 4000개
 *
 * 컴파일: gcc -O2 -Wall -o graph_compare graph_compare.c
 * 실행  : ./graph_compare [seed]      (기본 seed = 12345)
 *
 * 비교 횟수 기준:
 *   "연산을 위해 확인한 자료구조의 원소 수"
 *   - 인접 행렬  : matrix[a][b]를 한 번 확인할 때마다 1회
 *   - 인접 리스트: 리스트 노드를 한 개 확인할 때마다 1회
 *   리스트에 간선을 삽입할 때는 리스트의 맨 앞에 넣으므로, 중복 간선 확인 외에
 *   추가 비교가 없다. 삭제할 때는 양쪽 정점의 리스트에서 모두 노드를 찾아야 한다.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define V           100
#define SPARSE_E    100
#define DENSE_E     4000
#define TRIALS      1000   /* 평균을 구하기 위한 무작위 반복 횟수 */

typedef struct { int u, v; } Edge;

/* ---------- 인접 행렬 ---------- */
typedef struct { unsigned char m[V][V]; } AdjMatrix;

/* ---------- 인접 리스트 (연결 리스트) ---------- */
typedef struct Node {
    int vertex;
    struct Node *next;
} Node;

typedef struct {
    Node *head[V];
    int   nodeCount;          /* 할당된 간선 노드 수 (메모리 계산용) */
} AdjList;

static int rnd(int n) { return rand() % n; }

/* ================= 인접 행렬 ================= */
static int m_has(const AdjMatrix *g, int u, int v, long long *c) {
    (*c)++;
    return g->m[u][v];
}
static int m_insert(AdjMatrix *g, int u, int v, long long *c) {
    if (m_has(g, u, v, c)) return 0;       /* 이미 존재하는지 확인: 비교 1회 */
    g->m[u][v] = 1;
    g->m[v][u] = 1;
    return 1;
}
static int m_delete(AdjMatrix *g, int u, int v, long long *c) {
    if (!m_has(g, u, v, c)) return 0;      /* 존재 여부 확인: 비교 1회 */
    g->m[u][v] = 0;
    g->m[v][u] = 0;
    return 1;
}
static int m_adjacent(const AdjMatrix *g, int u, long long *c, int *out) {
    int n = 0;
    for (int j = 0; j < V; j++) {          /* 항상 100칸을 모두 확인 */
        (*c)++;
        if (g->m[u][j]) out[n++] = j;
    }
    return n;
}

/* ================= 인접 리스트 ================= */
static int l_has(const AdjList *g, int u, int v, long long *c) {
    for (const Node *p = g->head[u]; p; p = p->next) {
        (*c)++;
        if (p->vertex == v) return 1;
    }
    return 0;
}
static void push_front(AdjList *g, int u, int v) {
    Node *n = (Node *)malloc(sizeof(Node));
    n->vertex = v;
    n->next = g->head[u];
    g->head[u] = n;
    g->nodeCount++;
}
static int l_insert(AdjList *g, int u, int v, long long *c) {
    if (l_has(g, u, v, c)) return 0;       /* u의 리스트에서 중복 간선 확인 */
    push_front(g, u, v);                   /* 맨 앞에 삽입: 추가 비교 없음 */
    push_front(g, v, u);
    return 1;
}
/* u의 리스트에서 v를 삭제하며 확인한 노드 수를 센다. 삭제했으면 1 반환 */
static int remove_from(AdjList *g, int u, int v, long long *c) {
    Node **pp = &g->head[u];
    while (*pp) {
        (*c)++;
        if ((*pp)->vertex == v) {
            Node *del = *pp;
            *pp = del->next;
            free(del);
            g->nodeCount--;
            return 1;
        }
        pp = &(*pp)->next;
    }
    return 0;
}
static int l_delete(AdjList *g, int u, int v, long long *c) {
    if (!remove_from(g, u, v, c)) return 0;   /* u의 리스트에서 v를 찾아 삭제 */
    remove_from(g, v, u, c);                  /* v의 리스트에서 u를 찾아 삭제 */
    return 1;
}
static int l_adjacent(const AdjList *g, int u, long long *c, int *out) {
    int n = 0;
    for (const Node *p = g->head[u]; p; p = p->next) {
        (*c)++;
        out[n++] = p->vertex;
    }
    return n;
}
static void l_free(AdjList *g) {
    for (int i = 0; i < V; i++) {
        Node *p = g->head[i];
        while (p) { Node *t = p->next; free(p); p = t; }
        g->head[i] = NULL;
    }
    g->nodeCount = 0;
}

/* ================= 그래프 생성 ================= */
/* 중복과 자기 자신 연결이 없는 무작위 간선 ne개를 생성한다.
 * tmp는 난수 생성용 임시 배열이므로 그래프 메모리에 포함하지 않는다. */
static void gen_edges(Edge *e, int ne) {
    static unsigned char tmp[V][V];
    memset(tmp, 0, sizeof(tmp));
    int k = 0;
    while (k < ne) {
        int u = rnd(V), v = rnd(V);
        if (u == v || tmp[u][v]) continue;
        tmp[u][v] = tmp[v][u] = 1;
        e[k].u = u; e[k].v = v; k++;
    }
}

static long list_memory(const AdjList *g) {
    return (long)sizeof(AdjList) + (long)g->nodeCount * (long)sizeof(Node);
}

static int cmp_int(const void *a, const void *b) {
    return *(const int *)a - *(const int *)b;
}

/* 두 표현 방식이 실제로 같은 그래프인지 검증 */
static int verify_same(const AdjMatrix *m, const AdjList *l, int ne) {
    long long dummy = 0;
    int cnt = 0;
    for (int u = 0; u < V; u++) {
        if (m->m[u][u]) return 0;
        for (int v = 0; v < V; v++) {
            if (m->m[u][v] != m->m[v][u]) return 0;
            if (m->m[u][v] != l_has(l, u, v, &dummy)) return 0;
            if (m->m[u][v]) cnt++;
        }
    }
    return cnt == 2 * ne && l->nodeCount == 2 * ne;
}

static void print_list(const char *label, int u, int *a, int n) {
    printf("  %s 정점 %d (인접 정점 %d개): ", label, u, n);
    for (int i = 0; i < n; i++) printf("%d%s", a[i], i + 1 < n ? " " : "");
    printf("\n");
}

/* ================= 그래프 하나에 대한 실험 ================= */
static void run_graph(const char *name, int ne, int caseM, int caseL) {
    Edge *edges = (Edge *)malloc(sizeof(Edge) * ne);
    gen_edges(edges, ne);

    AdjMatrix *M = (AdjMatrix *)calloc(1, sizeof(AdjMatrix));
    AdjList   *L = (AdjList *)calloc(1, sizeof(AdjList));
    for (int i = 0; i < ne; i++) {
        M->m[edges[i].u][edges[i].v] = 1;
        M->m[edges[i].v][edges[i].u] = 1;
        push_front(L, edges[i].u, edges[i].v);
        push_front(L, edges[i].v, edges[i].u);
    }

    printf("==================== %s Graph (V=%d, E=%d) ====================\n", name, V, ne);
    printf("두 표현 방식이 같은 그래프인지 검증: %s\n\n",
           verify_same(M, L, ne) ? "OK" : "오류!");

    /* ----- 두 표현 방식에 공통으로 사용할 실험 대상 선택 ----- */
    long long dummy = 0;
    int iu, iv;                                  /* 삽입할 간선: 존재하지 않는 쌍 */
    do { iu = rnd(V); iv = rnd(V); } while (iu == iv || M->m[iu][iv]);
    Edge de = edges[rnd(ne)];                    /* 삭제할 간선: 존재하는 간선 */
    int cu, cv;                                  /* 연결 여부를 확인할 정점 쌍 */
    do { cu = rnd(V); cv = rnd(V); } while (cu == cv);
    int av = rnd(V);                             /* 인접 정점을 출력할 정점 */

    long long mi = 0, md = 0, mc = 0, ma = 0;
    long long li = 0, ld = 0, lc = 0, la = 0;
    int bufM[V], bufL[V];

    /* ----- 인접 행렬 ----- */
    long memM = (long)sizeof(AdjMatrix);
    m_insert(M, iu, iv, &mi);
    m_delete(M, de.u, de.v, &md);
    int hasM = m_has(M, cu, cv, &mc);
    int nM = m_adjacent(M, av, &ma, bufM);

    /* ----- 인접 리스트 (같은 순서로 연산) ----- */
    l_insert(L, iu, iv, &li);
    l_delete(L, de.u, de.v, &ld);
    int hasL = l_has(L, cu, cv, &lc);
    int nL = l_adjacent(L, av, &la, bufL);
    long memL = list_memory(L);

    printf("실험 대상 (두 방식에 동일하게 적용):\n");
    printf("  삽입 간선     : (%d, %d)  [존재하지 않음]\n", iu, iv);
    printf("  삭제 간선     : (%d, %d)  [존재함]\n", de.u, de.v);
    printf("  연결 확인     : (%d, %d)  -> %s\n", cu, cv, hasM ? "간선 있음" : "간선 없음");
    printf("  인접 정점 출력: 정점 %d\n\n", av);

    if (hasM != hasL) printf("!!! 두 방식의 연결 확인 결과가 다릅니다\n");

    /* 리스트 결과를 정렬한 사본으로 행렬 결과와 비교 */
    int sortedL[V];
    memcpy(sortedL, bufL, sizeof(int) * nL);
    qsort(sortedL, nL, sizeof(int), cmp_int);
    int same = (nM == nL) && memcmp(bufM, sortedL, sizeof(int) * nL) == 0;

    printf("Case %d : %s Graph - Adjacency Matrix\n", caseM, name);
    printf("Memory                         : %ld Bytes\n", memM);
    printf("Edge insertion comparisons     : %lld\n", mi);
    printf("Edge deletion comparisons      : %lld\n", md);
    printf("Connectivity check comparisons : %lld\n", mc);
    printf("Adjacent vertices comparisons  : %lld\n", ma);
    print_list("인접 정점:", av, bufM, nM);
    printf("\n");

    printf("Case %d : %s Graph - Adjacency List\n", caseL, name);
    printf("Memory                         : %ld Bytes\n", memL);
    printf("Edge insertion comparisons     : %lld\n", li);
    printf("Edge deletion comparisons      : %lld\n", ld);
    printf("Connectivity check comparisons : %lld\n", lc);
    printf("Adjacent vertices comparisons  : %lld\n", la);
    print_list("인접 정점:", av, bufL, nL);
    printf("\n두 방식의 인접 정점 집합 비교: %s\n\n", same ? "동일" : "다름!");

    /* 평균 측정을 위해 그래프를 원래 상태로 복구:
       삽입한 (iu,iv)를 삭제하고, 삭제했던 간선을 다시 삽입 */
    m_delete(M, iu, iv, &dummy);   l_delete(L, iu, iv, &dummy);
    m_insert(M, de.u, de.v, &dummy); l_insert(L, de.u, de.v, &dummy);

    /* ----- TRIALS회 무작위 실험의 평균 ----- */
    double sMi = 0, sMd = 0, sMc = 0, sMa = 0, sLi = 0, sLd = 0, sLc = 0, sLa = 0;
    for (int t = 0; t < TRIALS; t++) {
        long long a;
        int u, v;
        /* 삽입 */
        do { u = rnd(V); v = rnd(V); } while (u == v || M->m[u][v]);
        a = 0; m_insert(M, u, v, &a); sMi += a;
        a = 0; l_insert(L, u, v, &a); sLi += a;
        m_delete(M, u, v, &dummy); l_delete(L, u, v, &dummy);   /* 복구 */
        /* 삭제 */
        Edge e = edges[rnd(ne)];
        a = 0; m_delete(M, e.u, e.v, &a); sMd += a;
        a = 0; l_delete(L, e.u, e.v, &a); sLd += a;
        m_insert(M, e.u, e.v, &dummy); l_insert(L, e.u, e.v, &dummy); /* 복구 */
        /* 연결 확인 */
        do { u = rnd(V); v = rnd(V); } while (u == v);
        a = 0; m_has(M, u, v, &a); sMc += a;
        a = 0; l_has(L, u, v, &a); sLc += a;
        /* 인접 정점 */
        u = rnd(V);
        a = 0; m_adjacent(M, u, &a, bufM); sMa += a;
        a = 0; l_adjacent(L, u, &a, bufL); sLa += a;
    }
    printf("--- 무작위 %d회 실험 평균 (%s) ---\n", TRIALS, name);
    printf("%-26s %14s %14s\n", "연산", "인접 행렬", "인접 리스트");
    printf("%-26s %14.2f %14.2f\n", "간선 삽입",        sMi / TRIALS, sLi / TRIALS);
    printf("%-26s %14.2f %14.2f\n", "간선 삭제",        sMd / TRIALS, sLd / TRIALS);
    printf("%-26s %14.2f %14.2f\n", "연결 여부 확인",   sMc / TRIALS, sLc / TRIALS);
    printf("%-26s %14.2f %14.2f\n", "인접 정점 출력",   sMa / TRIALS, sLa / TRIALS);
    printf("\n");

    l_free(L);
    free(L); free(M); free(edges);
}

int main(int argc, char **argv) {
    unsigned seed = (argc > 1) ? (unsigned)atoi(argv[1]) : 12345u;
    srand(seed);
    printf("Seed = %u | sizeof(Node)=%zu, sizeof(AdjList)=%zu, sizeof(AdjMatrix)=%zu\n\n",
           seed, sizeof(Node), sizeof(AdjList), sizeof(AdjMatrix));

    run_graph("Sparse", SPARSE_E, 1, 2);
    run_graph("Dense",  DENSE_E,  3, 4);
    return 0;
}
