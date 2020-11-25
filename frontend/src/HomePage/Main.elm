module HomePage.Main exposing (..)

import Browser.Navigation as Nav
import Html exposing (..)
import Html.Attributes exposing (..)
import Html.Events exposing (onClick)
import Routes


type alias Model =
    { navKey : Nav.Key }


type Msg
    = SearchRuns
    | SearchProjects


init : Nav.Key -> Model
init navKey =
    { navKey = navKey }


update : Msg -> Model -> ( Model, Cmd msg )
update msg model =
    case msg of
        SearchRuns ->
            ( model, Nav.pushUrl model.navKey Routes.searchRunsRoute )

        SearchProjects ->
            ( model, Nav.pushUrl model.navKey Routes.searchProjectsRoute )


view : Html Msg
view =
    div [ class "d-sm-flex justify-content-center" ]
        [ div
            [ onClick SearchProjects
            , class "card m-4"
            , attribute "style" "width: 18rem;"
            , style "cursor" "pointer"
            ]
            [ img [ alt "Card image cap", class "card-img-top", src "project_card.png" ]
                []
            , div [ class "card-body" ]
                [ h5 [ class "card-title" ]
                    [ text "Search SRA ", b [] [ text "Projects" ] ]
                , p [ class "card-text" ]
                    [ text "Search for projects containing multiple sequencing runs" ]
                ]
            ]
        , div
            [ onClick SearchRuns
            , class "card m-4"
            , attribute "style" "width: 18rem;"
            , style "cursor" "pointer"
            ]
            [ img [ alt "Card image cap", class "card-img-top", src "run_card.png" ]
                []
            , div [ class "card-body " ]
                [ h5 [ class "card-title" ]
                    [ text "Search SRA ", b [] [ text "Runs" ] ]
                , p [ class "card-text" ]
                    [ text "Search for individual sequencing runs within a project" ]
                ]
            ]
        ]
