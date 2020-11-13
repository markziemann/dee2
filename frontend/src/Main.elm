module Main exposing (..)

import Browser exposing (Document)
import Browser.Navigation as Nav
import Elastic exposing (parse, serialize)
import Helpers exposing (decodeSearchResults)
import Html exposing (..)
import Html.Attributes exposing (..)
import Http
import Info exposing (introduction)
import Maybe.Extra
import Nav exposing (navbar)
import ResultsPage.Main as RPMain
import ResultsPage.Types
import ResultsPage.Views exposing (viewSearchResults)
import Routes
import SearchPage.Helpers exposing (getSearchMode, sameModeAndString, withPagination)
import SearchPage.Main as SPMain
import SearchPage.Types exposing (SearchMode(..), SearchParameters(..))
import SearchPage.Views exposing (viewLargeSearchBar, viewSearchButton, viewSearchModeSelector)
import SharedTypes exposing (PaginationOffset, RemoteData(..), toWebData)
import Types exposing (..)
import Url


init : () -> Url.Url -> Nav.Key -> ( Model, Cmd Msg )
init flags url navKey =
    let
        route =
            Routes.determinePage url

        model =
            { navKey = navKey
            , url = url
            , searchPage = SPMain.init navKey
            , resultsPage = RPMain.init navKey NotAsked Nothing
            , route = route
            }
    in
    case model.route of
        Routes.ResultsRoute searchUrlParameters ->
            ( model, search searchUrlParameters )

        _ ->
            ( model, Cmd.none )



---- UPDATE ----


search : SearchParameters -> Cmd Msg
search (SearchParameters searchMode searchString paginationOffset) =
    let
        parsedSearchString =
            case searchMode of
                Strict ->
                    Result.withDefault (Elastic.Word searchString) (parse searchString)
                        |> serialize

                Fuzzy ->
                    Elastic.Word searchString
                        |> serialize

        parsedParameters =
            SearchParameters searchMode parsedSearchString paginationOffset
    in
    Http.get
        { url = Routes.searchResultsRoute parsedParameters
        , expect =
            Http.expectJson
                (GotHttpSearchResponse parsedParameters << toWebData)
                (decodeSearchResults paginationOffset)
        }


update : Msg -> Model -> ( Model, Cmd Msg )
update msg model =
    let
        fromSearchPage =
            \( mdl, cmd ) ->
                ( { model | searchPage = mdl }, Cmd.map GotSearchPageMsg cmd )

        fromResultsPage =
            \( mdl, cmd ) ->
                ( { model | resultsPage = mdl }, Cmd.map GotResultsPageMsg cmd )
    in
    case msg of
        GotSearchPageMsg message ->
            SPMain.update message model.searchPage
                |> fromSearchPage

        GotResultsPageMsg message ->
            RPMain.update message model.resultsPage
                |> fromResultsPage

        LinkClicked urlRequest ->
            case urlRequest of
                Browser.Internal url ->
                    Debug.log "Second" ( model, Nav.pushUrl model.navKey (Url.toString url) )

                Browser.External href ->
                    ( model, Nav.load href )

        GotHttpSearchResponse searchParameters webData ->
            ( model, Cmd.none )

        UrlChanged url ->
            let
                route =
                    Routes.determinePage url

                default =
                    \cmd ->
                        ( { model | url = url, route = route }, cmd )
            in
            case route of
                Routes.ResultsRoute searchParameters ->
                    default <| search searchParameters

                _ ->
                    default Cmd.none



---- VIEW ----


pageLayout : List (Html Msg) -> List (Html Msg)
pageLayout content =
    [ navbar
    , div [ class "container my-5 mx-auto" ] content
    , introduction
    ]


pageView : Model -> List (Html Msg)
pageView model =
    let
        fromSearchPage =
            List.map (Html.map GotSearchPageMsg)

        fromResultsPage =
            List.map (Html.map GotResultsPageMsg)
    in
    pageLayout <|
        case model.route of
            Routes.HomeRoute ->
                fromSearchPage
                    [ viewLargeSearchBar model.searchPage
                    , viewSearchModeSelector <| getSearchMode model.searchPage.searchParameters
                    , viewSearchButton model.searchPage
                    ]

            Routes.ResultsRoute (SearchParameters _ _ paginationOffset) ->
                fromResultsPage
                    (viewSearchResults model.resultsPage paginationOffset)

            Routes.Unknown ->
                [ text "Hmm... I don't recognise that url." ]


view : Model -> Document Msg
view model =
    { title = "Digital Expression Explorer 2"
    , body = pageView model
    }


subscriptions : Model -> Sub Msg
subscriptions model =
    case model.route of
        Routes.HomeRoute ->
            Sub.batch
                [ Sub.map GotSearchPageMsg <| SPMain.subscriptions model.searchPage
                ]

        Routes.ResultsRoute _ ->
            Sub.none

        Routes.Unknown ->
            subscriptions { model | route = Routes.HomeRoute }



---- PROGRAM ----


main : Program () Model Msg
main =
    Browser.application
        { view = view
        , init = init
        , update = update
        , subscriptions = subscriptions
        , onUrlRequest = LinkClicked
        , onUrlChange = UrlChanged
        }
